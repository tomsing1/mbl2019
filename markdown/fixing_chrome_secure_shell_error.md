# Introduction

Sometimes the Chrome Secure Shell shows an error like this:

```
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ WARNING: REMOTE HOST IDENTIFICATION HAS CHANGED! @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
IT IS POSSIBLE THAT SOMEONE IS DOING SOMETHING NASTY!
Someone could be eavesdropping on you right now (man-in-the-middle attack)!
It is also possible that a host key has just been changed.
The fingerprint for the ECDSA key sent by the remote host is
aa:bb:cc:dd:ee:ff:gg:hh:ii:jj:kk:ll:mm:nn:oo:pp.
[...]
```

To fix it:

1. Try to establish a connection and when you get the error press the `C` key to get back the input form.
2. While the form is active, open the Javascript console by pressing the following keys:

- Windows: CTRL +Shift + j
- Mac: OPTION + COMMAND + j
3. Paste the following code into the Javascript console and press enter:

```
term_.command.removeAllKnownHosts()
```

