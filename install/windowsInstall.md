# Windows Setup

For those running Windows, we could probably perform about half of today's session on your computer.
However, as many tools rely on the fundamental architecture of OSX & Linux systems, not all tools are available.
The best solution for today's session is to use a Virtual Machine (VM) which we have configured for you.
This is essentially a computer running somewhere else (i.e. the "*cloud*") which we can log onto remotely.

To setup a connection to this VM we first need to install the client software X2GO, which can be obtained [from here](https://code.x2go.org/releases/binary-win32/x2goclient/releases/4.1.0.0-2017.03.11/x2goclient-4.1.0.0-2017.03.11-setup.exe).
Once you've installed this software, **ask the instructors for an IP address and password.**

To connect to your VM, please follow these instructions carefully.
First, we need to create a session with the basic parameters, then we'll actually log in.

1. Open X2GO
2. Enter *Intro-NGS* as the Session Name
3. Enter your *IP address* where it says *Host*
4. Enter the word `ubuntu` as the login. **This must be all lower-case**
5. Select XFCE from the drop-down menu under Session Type
6. Click OK

Now we have created the session, it will appear in X2Go on the right panel.
To log onto the VM, we simply click on the session, and enter the password you were given.
Click OK if you receive a message about a security key.
*If this process fails, please place a post-it note on your monitor.*

**We advise maximising your X2Go window to replicate sitting at the VM as if it is your local machine.**

Once you've logged in you should see a desktop like the one below.
To open the terminal, click the *terminal icon* at the bottom of the screen.

![](images/VM_Desktop.png)
