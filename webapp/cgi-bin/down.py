#!/usr/bin/python3
# encoding : utf-8
 
import flup.server.fcgi as flups #fast cgi
from flup.server.fcgi import WSGIServer
from email import encoders
from email.header import Header
from email.mime.text import MIMEText
from email.utils import parseaddr, formataddr
import smtplib
from urllib.parse import unquote
import json
import sys
import os
import requests
def get_environ(environ):
    # rquest_method = environ["REQUEST_METHOD"]
    # str = "rquest_method:" + rquest_method + "\r\n"
    # query_string = environ["QUERY_STRING"]
    # str += ",query_string:" + query_string + "\r\n"
    # script_filename = environ["SCRIPT_FILENAME"]
    # str += ",script_filename:" + script_filename + "\r\n"
    # script_name = environ["SCRIPT_NAME"]
    # str += ",script_name:" + script_name + "\r\n"
    # rquest_uri = environ["REQUEST_URI"]
    # str += ", rquest_uri:" + rquest_uri + "\r\n"
    # remote_addr = environ["REMOTE_ADDR"]
    # str += ",remote_addr:" + remote_addr + "\r\n"
    # remote_port = environ["REMOTE_PORT"]
    # str += ",remote_port:" + remote_port + "\r\n"
    
    data = environ["wsgi.input"].read()
    str1 = bytes.decode(data)
    return str1 

def work(gene):
    stat=os.system("""
    cd /var/www/hap.bioinf.club/webapp/gene_data;
    tar cvhzf %s.tar.gz %s*.fasta %s*.csv 
    """%(gene,gene,gene))

    down_url=''
    if not stat:
        down_url="http://hap.bioinf.club/gene_data/%s.tar.gz"%(gene)
    
    return tuple([stat,down_url])




def myapp(environ, start_response):
    content = get_environ(environ)
    #name=xiao&email=fengcong97%40qq.com&subject=test&message=1111111111
    split_l = content.split("&")
    
    gene= unquote(split_l[0].split("name=")[1],'utf-8')
    act= unquote(split_l[1].split("act=")[1],'utf-8')
    stat=0
    down_url = ''

    if act == "down":
        stat,down_url = work(gene)
    os.system("echo stat=%d,down_url=%s > /var/www/hap.bioinf.club/webapp/gene_data/test.txt"%(stat,down_url))
    #os.system("echo 11 > /var/www/hap.bioinf.club/webapp/gene_data/test.out")   
    start_response('200 OK', [('Content-Type', 'text/plain')])
    
    if not stat :
        return down_url
    else:
        return "notOK"

    
    

if __name__  == '__main__':
    #WSGIServer(myapp,bindAddress=('127.0.0.1',8008)).run()
    flups.WSGIServer(myapp).run() #fastcgi
