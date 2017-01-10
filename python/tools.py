
# -*- coding: utf-8 -*-

import subprocess


def call_command(*args, **kwargs):
	popen_kwargs = {
		"stdin" : subprocess.PIPE,
		"stdout" : subprocess.PIPE,
		"stderr" : subprocess.PIPE,
		"shell" : True,
	}
	popen_kwargs.update(kwargs)
	process = subprocess.Popen(args, **popen_kwargs)
	return process.communicate()

def convert_certificate_key_rsa(in_key, out_key):
	call_command("openssl rsa -in {in_key} -out {out_key}".format(in_key=in_key, out_key=out_key))

