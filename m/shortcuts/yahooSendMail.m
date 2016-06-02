function yahooSendMail(myMail,MyPassword,recipients,subject,message,attachments)

if nargin<6
    
    attachments={};
    
end

setpref('Internet','SMTP_Server','smtp.mail.yahoo.com');

setpref('Internet','E_mail',myMail);

setpref('Internet','SMTP_Username',myMail);

setpref('Internet','SMTP_Password',MyPassword);

props = java.lang.System.getProperties;

props.setProperty('mail.smtp.auth','true');

props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');

props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(recipients,subject,message,attachments);

end

