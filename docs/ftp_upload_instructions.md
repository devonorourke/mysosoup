installed FTP with brew.
went to SRA upload site and selected `FTP upload` which had pre-generated tmp username/password info
then used Terminal to connect to FTP site, starting at directory where files were downloaded

```
ftp -i ftp-private.ncbi.nlm.nih.gov   ## logs into FTP site
subftp    ## username
{enter password}  ## see drop down menu on Submission site for details
cd uploads/devon.orourke_gmail.com_SFmljtHQ   ## go to subfolder
mkdir mysosoup    ## make subfolder within this subfolder
cd mysosoup     ## move into subdirectory
passive on  ## avoids firewall issues
mput *.gz    ## upload all the .gz files (don't just use 'put')
```
This will work; just takes forever.
