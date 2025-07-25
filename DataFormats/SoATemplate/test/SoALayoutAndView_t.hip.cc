#include <cassert>
#include <cstdlib>
#include <memory>

#include <hip/hip_runtime.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "HeterogeneousCore/ROCmUtilities/interface/hipCheck.h"
#include "HeterogeneousCore/ROCmUtilities/interface/requireDevices.h"

// Test SoA stores and view.
// Use cases
// Multiple stores in a buffer
// Scalars, Columns of scalars and of Eigen vectors
// View to each of them, from one and multiple stores.

GENERATE_SOA_LAYOUT(SoAHostDeviceLayoutTemplate,
                    /*SoAHostDeviceViewTemplate,*/
                    // predefined static scalars
                    // size_t size;
                    // size_t alignment;

                    // columns: one value per element
                    SOA_COLUMN(double, x),
                    SOA_COLUMN(double, y),
                    SOA_COLUMN(double, z),
                    SOA_EIGEN_COLUMN(Eigen::Vector3d, a),
                    SOA_EIGEN_COLUMN(Eigen::Vector3d, b),
                    SOA_EIGEN_COLUMN(Eigen::Vector3d, r),
                    // scalars: one value for the whole structure
                    SOA_SCALAR(const char*, description),
                    SOA_SCALAR(uint32_t, someNumber))

using SoAHostDeviceLayout = SoAHostDeviceLayoutTemplate<>;
using SoAHostDeviceView = SoAHostDeviceLayout::View;
using SoAHostDeviceRangeCheckingView =
    SoAHostDeviceLayout::ViewTemplate<cms::soa::RestrictQualify::enabled, cms::soa::RangeChecking::enabled>;
using SoAHostDeviceConstView = SoAHostDeviceLayout::ConstView;

GENERATE_SOA_LAYOUT(SoADeviceOnlyLayoutTemplate,
                    /*SoADeviceOnlyViewTemplate,*/
                    SOA_COLUMN(uint16_t, color),
                    SOA_COLUMN(double, value),
                    SOA_COLUMN(double*, py),
                    SOA_COLUMN(uint32_t, count),
                    SOA_COLUMN(uint32_t, anotherCount))

using SoADeviceOnlyLayout = SoADeviceOnlyLayoutTemplate<>;
using SoADeviceOnlyView = SoADeviceOnlyLayout::View;

GENERATE_SOA_LAYOUT(SoAFullDeviceLayoutTemplate,
                    SOA_COLUMN(double, x),
                    SOA_COLUMN(double, y),
                    SOA_COLUMN(double, z),
                    SOA_COLUMN(uint16_t, color),
                    SOA_COLUMN(double, value),
                    SOA_COLUMN(double*, py),
                    SOA_COLUMN(uint32_t, count),
                    SOA_COLUMN(uint32_t, anotherCount),
                    SOA_SCALAR(const char*, description),
                    SOA_SCALAR(uint32_t, someNumber))

using SoAFullDeviceLayout =
    SoAFullDeviceLayoutTemplate<cms::soa::CacheLineSize::NvidiaGPU, cms::soa::AlignmentEnforcement::enforced>;
using SoAFullDeviceView = SoAFullDeviceLayout::View;
using SoAFullDeviceConstView = SoAFullDeviceLayout::ConstView;

// These SoAs validate that the generating macros do not get confused in the special case where there are
// no columns and only scalar elements in the SoA.
GENERATE_SOA_LAYOUT(TestSoALayoutNoColumn, SOA_SCALAR(double, r))
GENERATE_SOA_LAYOUT(TestSoALayoutNoColumn2, SOA_SCALAR(double, r), SOA_SCALAR(double, r2))

// Eigen cross product kernel (on store)
__global__ void crossProduct(SoAHostDeviceView soa, const int numElements) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= numElements)
    return;
  auto si = soa[i];
  si.r() = si.a().cross(si.b());
}

// Device-only producer kernel
__global__ void producerKernel(SoAFullDeviceView soa, const int numElements) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= numElements)
    return;
  auto si = soa[i];
  si.color() &= 0x55 << i % (sizeof(si.color()) - sizeof(char));
  si.value() = sqrt(si.x() * si.x() + si.y() * si.y() + si.z() * si.z());
}

// Device-only consumer with result in host-device area
__global__ void consumerKernel(SoAFullDeviceView soa, const int numElements) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= numElements)
    return;
  auto si = soa[i];
  si.x() = si.color() * si.value();
}

// Get a view like the default, except for range checking
using RangeCheckingHostDeviceView =
    SoAHostDeviceLayout::ViewTemplate<SoAHostDeviceView::restrictQualify, cms::soa::RangeChecking::enabled>;

// We expect to just run one thread.
__global__ void rangeCheckKernel(RangeCheckingHostDeviceView soa) {
  printf("About to fail range-check (operator[]) in ROCm thread: %d\n", (int)threadIdx.x);
  [[maybe_unused]] auto si = soa[soa.metadata().size()];
  printf("Fail: range-check failure should have stopped the kernel.\n");
}

int main(void) {
  cms::rocmtest::requireDevices();

  hipStream_t stream;
  hipCheck(hipStreamCreateWithFlags(&stream, hipStreamNonBlocking));

  // Non-aligned number of elements to check alignment features.
  constexpr unsigned int numElements = 65537;

  // Allocate buffer and store on host
  size_t hostDeviceSize = SoAHostDeviceLayout::computeDataSize(numElements);
  std::byte* h_buf = nullptr;
  hipCheck(hipHostMalloc((void**)&h_buf, hostDeviceSize));
  SoAHostDeviceLayout h_soahdLayout(h_buf, numElements);
  SoAHostDeviceView h_soahd(h_soahdLayout);

  // Validation of range checking variants initialization
  SoAHostDeviceRangeCheckingView h_soahdrc(h_soahdLayout);
  [[maybe_unused]] SoAHostDeviceRangeCheckingView h_soahdrc2 = h_soahdLayout;
  [[maybe_unused]] SoAHostDeviceRangeCheckingView h_soahdrc3{h_soahd};
  [[maybe_unused]] SoAHostDeviceRangeCheckingView h_soahdrc4 = h_soahd;
  SoAHostDeviceConstView h_soahd_c(h_soahdLayout);

  // Alocate buffer, stores and views on the device (single, shared buffer).
  size_t deviceOnlySize = SoADeviceOnlyLayout::computeDataSize(numElements);
  std::byte* d_buf = nullptr;
  hipCheck(hipHostMalloc((void**)&d_buf, hostDeviceSize + deviceOnlySize));
  SoAHostDeviceLayout d_soahdLayout(d_buf, numElements);
  SoADeviceOnlyLayout d_soadoLayout(d_soahdLayout.metadata().nextByte(), numElements);
  SoAHostDeviceView d_soahdView(d_soahdLayout);
  SoADeviceOnlyView d_soadoView(d_soadoLayout);
  const auto d_soahdRecords = d_soahdView.records();
  const auto d_soadoRecords = d_soadoView.records();
  SoAFullDeviceView d_soaFullView(d_soahdRecords.x(),
                                  d_soahdRecords.y(),
                                  d_soahdRecords.z(),
                                  d_soadoRecords.color(),
                                  d_soadoRecords.value(),
                                  d_soadoRecords.py(),
                                  d_soadoRecords.count(),
                                  d_soadoRecords.anotherCount(),
                                  d_soahdRecords.description(),
                                  d_soahdRecords.someNumber());

  // Assert column alignments
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_x()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_y()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_z()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_a()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_b()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_r()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_description()) % decltype(h_soahd)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(h_soahd.metadata().addressOf_someNumber()) % decltype(h_soahd)::alignment);

  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_x()) % decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_y()) % decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_z()) % decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_a()) % decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_b()) % decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_r()) % decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_description()) %
                  decltype(d_soahdLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soahdLayout.metadata().addressOf_someNumber()) %
                  decltype(d_soahdLayout)::alignment);

  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soadoLayout.metadata().addressOf_color()) % decltype(d_soadoLayout)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soadoLayout.metadata().addressOf_value()) % decltype(d_soadoLayout)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soadoLayout.metadata().addressOf_py()) % decltype(d_soadoLayout)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soadoLayout.metadata().addressOf_count()) % decltype(d_soadoLayout)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soadoLayout.metadata().addressOf_anotherCount()) %
                  decltype(d_soadoLayout)::alignment);

  // Views should get the same alignment as the stores they refer to
  assert(0 == reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_x()) % decltype(d_soaFullView)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_y()) % decltype(d_soaFullView)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_z()) % decltype(d_soaFullView)::alignment);
  // Limitation of views: we have to get scalar member addresses via metadata.
  assert(0 == reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_description()) %
                  decltype(d_soaFullView)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_someNumber()) %
                  decltype(d_soaFullView)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_color()) % decltype(d_soaFullView)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_value()) % decltype(d_soaFullView)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_py()) % decltype(d_soaFullView)::alignment);
  assert(0 ==
         reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_count()) % decltype(d_soaFullView)::alignment);
  assert(0 == reinterpret_cast<uintptr_t>(d_soaFullView.metadata().addressOf_anotherCount()) %
                  decltype(d_soaFullView)::alignment);

  // Initialize and fill the host buffer
  std::memset(h_soahdLayout.metadata().data(), 0, hostDeviceSize);
  for (size_t i = 0; i < numElements; ++i) {
    auto si = h_soahd[i];
    // Tuple assignment...
    // elements are: x, y, z, a, b, r
    auto v1 = 1.0 * i + 1.0;
    auto v2 = 2.0 * i;
    auto v3 = 3.0 * i - 1.0;
    if (i % 2) {
      si = {v1, v2, v3, {v1, v2, v3}, {v3, v2, v1}, {0, 0, 0}};
    } else {
      si.x() = si.a()(0) = si.b()(2) = v1;
      si.y() = si.a()(1) = si.b()(1) = v2;
      si.z() = si.a()(2) = si.b()(0) = v3;
    }
  }
  auto& sn = h_soahd.someNumber();
  sn = numElements + 2;

  // Push to device
  hipCheck(hipMemcpyAsync(d_buf, h_buf, hostDeviceSize, hipMemcpyDefault, stream));

  // Process on device
  crossProduct<<<(numElements + 255) / 256, 256, 0, stream>>>(d_soahdView, numElements);

  // Paint the device only with 0xFF initially
  hipCheck(hipMemsetAsync(d_soadoLayout.metadata().data(), 0xFF, d_soadoLayout.metadata().byteSize(), stream));

  // Produce to the device only area
  producerKernel<<<(numElements + 255) / 256, 256, 0, stream>>>(d_soaFullView, numElements);

  // Consume the device only area and generate a result on the host-device area
  consumerKernel<<<(numElements + 255) / 256, 256, 0, stream>>>(d_soaFullView, numElements);

  // Get result back
  hipCheck(hipMemcpyAsync(h_buf, d_buf, hostDeviceSize, hipMemcpyDefault, stream));

  // Wait and validate.
  hipCheck(hipStreamSynchronize(stream));
  for (size_t i = 0; i < numElements; ++i) {
    auto si = h_soahd_c[i];
    assert(si.r() == si.a().cross(si.b()));
    double initialX = 1.0 * i + 1.0;
    double initialY = 2.0 * i;
    double initialZ = 3.0 * i - 1.0;
    uint16_t expectedColor = 0x55 << i % (sizeof(uint16_t) - sizeof(char));
    double expectedX = expectedColor * sqrt(initialX * initialX + initialY * initialY + initialZ * initialZ);
    if (abs(si.x() - expectedX) / expectedX >= 2 * std::numeric_limits<double>::epsilon()) {
      std::cout << "X failed: for i=" << i << std::endl
                << "initialX=" << initialX << " initialY=" << initialY << " initialZ=" << initialZ << std::endl
                << "expectedX=" << expectedX << std::endl
                << "resultX=" << si.x() << " resultY=" << si.y() << " resultZ=" << si.z() << std::endl
                << "relativeDiff=" << abs(si.x() - expectedX) / expectedX
                << " epsilon=" << std::numeric_limits<double>::epsilon() << std::endl;
      assert(false);
    }
  }

  {
    // Get a view like the default, except for range checking (direct initialization from layout)
    SoAHostDeviceRangeCheckingView soa1viewRangeChecking(h_soahdLayout);
    try {
      [[maybe_unused]] auto si = soa1viewRangeChecking[soa1viewRangeChecking.metadata().size()];
      std::cout << "Fail: expected range-check exception (view-level index access) not caught on the host (overflow)."
                << std::endl;
      assert(false);
    } catch (const std::out_of_range&) {
    }
    try {
      [[maybe_unused]] auto si = soa1viewRangeChecking[-1];
      std::cout << "Fail: expected range-check exception (view-level index access) not caught on the host (underflow)."
                << std::endl;
      assert(false);
    } catch (const std::out_of_range&) {
    }
    [[maybe_unused]] auto si = soa1viewRangeChecking[soa1viewRangeChecking.metadata().size() - 1];
    [[maybe_unused]] auto si2 = soa1viewRangeChecking[0];
    std::cout << "Pass: expected range-check exceptions (view-level index access) successfully caught on the host "
                 "(layout initialization)."
              << std::endl;
  }

  {
    // Validation of view initialized range checking view initialization
    try {
      [[maybe_unused]] auto si = h_soahdrc3[h_soahdrc3.metadata().size()];
      std::cout << "Fail: expected range-check exception (view-level index access) not caught on the host (overflow)."
                << std::endl;
      assert(false);
    } catch (const std::out_of_range&) {
    }
    try {
      [[maybe_unused]] auto si = h_soahdrc3[-1];
      std::cout << "Fail: expected range-check exception (view-level index access) not caught on the host (underflow)."
                << std::endl;
      assert(false);
    } catch (const std::out_of_range&) {
    }
    [[maybe_unused]] auto si = h_soahdrc3[h_soahdrc3.metadata().size() - 1];
    [[maybe_unused]] auto si2 = h_soahdrc3[0];
    std::cout << "Pass: expected range-check exceptions (view-level index access) successfully caught on the host "
                 "(view initialization)."
              << std::endl;
  }

  // Validation of range checking in a kernel
  // Disable this test until ROCm provides a non-fatal way to assert in device code
#if 0
  // Get a view like the default one, except for range checking
  RangeCheckingHostDeviceView soa1viewRangeChecking(d_soahdLayout);

  // This should throw an exception in the kernel
  rangeCheckKernel<<<1, 1, 0, stream>>>(soa1viewRangeChecking);

  // Wait and confirm that the ROCm kernel failed
  try {
    hipCheck(hipStreamSynchronize(stream));
    std::cout << "Fail: expected range-check exception not caught while executing the kernel." << std::endl;
    assert(false);
  } catch (const std::runtime_error&) {
    std::cout << "Pass: expected range-check exception caught while executing the kernel." << std::endl;
  }
#endif

  std::cout << "OK" << std::endl;
}
