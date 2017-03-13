
# -*- coding: utf-8 -*-

import multiprocessing
import subprocess
import time


def call_command(*args, **kwargs):
	print " ".join(args)
	popen_kwargs = {
		"stdin" : subprocess.PIPE,
		"stdout" : subprocess.PIPE,
		"stderr" : subprocess.PIPE,
		"shell" : True,
	}
	popen_kwargs.update(kwargs)
	process = subprocess.Popen(args, **popen_kwargs)
	return process.communicate()

def parallelize(function, arguments_list, n_processes=1):
	if n_processes <= 1:
		results = []
		for arguments in arguments_list:
			results.append(function(arguments))
		return results
	else:
		pool = multiprocessing.Pool(processes=max(1, min(n_processes, len(arguments_list))))
		results = pool.map_async(function, arguments_list)
		n_tasks = len(arguments_list)
		left = n_tasks
		iterator = iter(range(n_tasks))
		iterator.next()
		while (True):
			if (results.ready()): break
			remaining = results._number_left
			if remaining < left:
				for i in range(left-remaining):
					iterator.next()
				left = remaining
			time.sleep(1.0)
		returnvalue = results.get(9999999)
		pool.close() # necessary to actually terminate the processes
		pool.join()  # without these two lines, they happen to live until the whole program terminates
		return returnvalue

def convert_certificate_key_rsa(in_key, out_key):
	call_command("openssl rsa -in {in_key} -out {out_key}".format(in_key=in_key, out_key=out_key))

