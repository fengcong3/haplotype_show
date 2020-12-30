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



def myapp(environ, start_response):
    content = get_environ(environ)
    #name=xiao&email=fengcong97%40qq.com&subject=test&message=1111111111
    split_l = content.split("&")

    gene= unquote(split_l[0].split("gene=")[1],'utf-8')

    ret_json = {}
    try:
        with open("../gene_data/%s.json",'r') as load_f:
            ret_json = json.load(load_f)   
    except Exception:
        ret_json = {}
     

    start_response('200 OK', [('Content-Type', 'text/plain')])

    return ret_json

    
    

if __name__  == '__main__':
    #WSGIServer(myapp,bindAddress=('127.0.0.1',8008)).run()
    flups.WSGIServer(myapp).run() #fastcgi
