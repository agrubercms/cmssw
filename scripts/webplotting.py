#!/usr/bin/env python

import argparse

import Validation.tau_validation_tools.webplotting as webplotting


def main():
	parser = argparse.ArgumentParser(description="Create web gallery for a set of plots.")

	parser.add_argument("input_dir", help="Input directory containing the plots")
	parser.add_argument("-r", "--recursive", default=False, action="store_true",
	                    help="Recursively search for plots in the input directory.")
	parser.add_argument("--file-types", default=["png", "jpg", "pdf", "svg", "eps"],
	                    help="File types to be considered as plots. [Default: %(default)s]")

	args = parser.parse_args()
	webplotting.webplotting(
			input_dir=args.input_dir,
			recursive=args.recursive,
			file_types=args.file_types
	)
	

if __name__ == "__main__":
	main()
