INSTALLATION

1. Create a directory which will house the components from the download
   package or cvs, eg:
     % mkdir $HOME/components

2. Set your COMPONENTS_HOME environment variable to that location (you should
   add this to your .login or .post.cshrc file):
     % setenv COMPONENTS_HOME $HOME/components

3. Choose or create a bin/ directory where the executables will reside.  This 
   directory should be considered volatile, ie, only housing things that
   are installed and can be deleted by an uninstall.  It should not hold
   hand-edited files. You can place it in your home or a public place, eg:
     % mkdir $HOME/bin

4. Add this directory to the front of your PATH (you should
   add this to your .login or .post.cshrc file):
     % setenv PATH $HOME/bin:$PATH

5. Choose or create a lib/perl/ directory where the perl libraries will 
   reside. This directory should be considered volatile, ie, only housing 
   things that are installed and can be deleted by an uninstall.  It should 
   not hold hand-edited files.  You can place it in your home or a public 
   place, eg:
     % mkdir -p $HOME/lib/perl

6. Add that directory to the front of your PERL5LIB (you should
   add this to your .login or .post.cshrc file):
     % setenv PERL5LIB $HOME/lib/perl:$PERL5LIB

7. Create the test/ directory where the test data will reside.  This 
   directory should be considered volatile, ie, only housing things that
   are installed and can be deleted by an uninstall.  It should not hold
   hand-edited files.  You can place it in your home or a public place. THE 
   DEFAULT LOCATION IS YOUR HOME.  If you place the test/ directory there, 
   you will need to do less reconfiguring to run the tests.
     % mkdir -p $HOME/test

If you have access to CBIL's cvs repository:

   1. checkout the DistribJobTasks install script
     % cd $COMPONENTS_HOME
     % cvs co DistribJobTasks/install

   2. Run the install script (example below assumes you placed bin/ and lib/
      in $HOME).
     % cd $COMPONENTS_HOME
     % DistribJobTasks/install $HOME/bin $HOME/lib/perl $HOME/test -checkoutfirst


If you are using a download package:

   1. Run the install script (example below assumes you placed bin/ and lib/
      in $HOME).
     % cd $COMPONENTS_HOME
     % tar -xf the_download_file
     % DistribJobTasks/install $HOME/bin $HOME/lib/perl $HOME/test



