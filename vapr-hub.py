#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 12:07:36 2017

@author: nima9589
"""

import paramiko

VHOST = 'xenon.colorado.edu'
VHKEY = paramiko.hostkeys.HostKeyEntry.from_line('|1|Ibuxn8d47IJEVVmr/iuKPef/cvs=|SsfkMWNAq9p+rRLJL4m0WvXJu7k= ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA2Mca2FB7zjravXTzPwK4Cdy9dzTnSbAVmAQXUzAZT8XHoPyviENPGQCJv1YiGX4eoohPqtSLLzD6Ktz/ONf9/y89cDGNtzAX/VnVX5ejQWbd6DMd+6PwazEEW9Ccf4tKgXb9tlllTqSQ5MJVBltmBtF9ULNZMMyHvIGslJFe05zYOBNwx3VPCVf9YL3xexUPP6jcYZMiZ+KTxMyyUWmkpwfmek596/Jier4J3Rh/czOQ90yd+7lBKrLHp3JtC6Tv7oVgazNRU+SdK4x/PeBVWs8MEmG1uAWwMly1dOd//IeAAmuca56hgHBWCTznerTOoIjbmOa1I1nAR5TXLbWWDQ==')
VUSER = 'vapr'
VPKEY = paramiko.pkey.PKey.from_private_key_file('/home/vapr/.ssh/id_vapr')
VCMD = 'vapr-server'

with paramiko.transport.Transport(VHOST) as transpo:
    #transpo.use_compression(True)
    transpo.connect(hostkey=VHKEY, username=VUSER, pkey=VPKEY)
    with transpo.open_session() as chan:
        chan.exec_command(VCMD)
        byt = chan.recv(50)
	if not byt:
            # channel is closed
        # do stuff with byt
        sendall(data) # blocks until all of data is sent
        # send(data) # just sends some, and returns number of bytes sent
