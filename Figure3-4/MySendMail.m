function MySendMail

MailAddress = 'shiymqh@gmail.com';
password = '602672qq';  
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
subject = 'From Yuanming MATLAB';
content = 'OK!';
sendmail('shiymqh@gmail.com',subject,content);
end

