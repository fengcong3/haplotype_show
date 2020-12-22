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

def mail(name,email,subject,message): 
    my_sender='g3rp@qq.com' # 发件人邮箱账号 
    my_pass = 'pketbcounxdbcbac' # 发件人邮箱密码 
    my_user='fengcong97@qq.com' # 收件人邮箱账号，我这边发送给自己 
    ret=True 
    try:
        txt=subject+"\r\n\r\n"+message+"\r\n\r\nfrom: "+name+"("+email+")"
        msg=MIMEText(txt,'plain','utf-8')       
        msg['From']=formataddr(["G3RP",my_sender]) 
        # 括号里的对应发件人邮箱昵称、发件人邮箱账号 
        msg['To']=formataddr(["fengcong",my_user]) 
        # 括号里的对应收件人邮箱昵称、收件人邮箱账号 
        msg['Subject']="G3RP MESSAGE"
        # # 邮件的主题，也可以说是标题 
        server=smtplib.SMTP_SSL("smtp.qq.com", 465) 
        # # 发件人邮箱中的SMTP服务器，端口是25 
        server.login(my_sender, my_pass) 
        # # 括号中对应的是发件人邮箱账号、邮箱密码 
        server.sendmail(my_sender,[my_user,],msg.as_string()) 
        # # 括号中对应的是发件人邮箱账号、收件人邮箱账号、发送邮件 
        server.quit() # 关闭连接 
    except Exception: 
	    # 如果 try 中的语句没有执行，则会执行下面的 ret=False 
        ret=False
    
    return ret

def myapp(environ, start_response):
    content = get_environ(environ)
    #name=xiao&email=fengcong97%40qq.com&subject=test&message=1111111111
    split_l = content.split("&")

    name= unquote(split_l[0].split("name=")[1],'utf-8')
    email= split_l[1].split("email=")[1].replace("%40","@")
    subject=unquote(split_l[2].split("subject=")[1].replace("+"," "),'utf-8')
    message = unquote(split_l[3].split("message=")[1].replace("+"," "),'utf-8')
    

    
     
    ret=mail(name,email,subject,message) 
    start_response('200 OK', [('Content-Type', 'text/plain')])

    if ret: 
        return ['Your message has been sent. Thank you!'] 
    else: 
        return ['Server error!']

    
    

if __name__  == '__main__':
    #WSGIServer(myapp,bindAddress=('127.0.0.1',8008)).run()
    flups.WSGIServer(myapp).run() #fastcgi
