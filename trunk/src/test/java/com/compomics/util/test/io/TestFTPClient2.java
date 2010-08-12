/*
 * Copyright (C) Lennart Martens
 * 
 * Contact: lennart.martens AT UGent.be (' AT ' to be replaced with '@')
 */

/*
 * Created by IntelliJ IDEA.
 * User: Lennart
 * Date: 25-nov-02
 * Time: 11:21:52
 */
package com.compomics.util.test.io;
import org.apache.log4j.Logger;

import junit.TestCaseLM;

import java.util.Properties;
import java.io.*;

import com.compomics.util.io.FTPClient2;
import junit.framework.*;

/*
 * CVS information:
 *
 * $Revision: 1.3 $
 * $Date: 2007/07/06 09:41:53 $
 */

/**
 * This class implements the test scenario for an FTPClient2.
 *
 * @author Lennart Martens
 * @see com.compomics.util.io.FTPClient
 */
public class TestFTPClient2 extends TestCaseLM {

    // Class specific log4j logger for TestFTPClient2 instances.
    Logger logger = Logger.getLogger(TestFTPClient2.class);

    private String iServer = null;
    private String iUser = null;
    private String iPassword = null;

    private boolean iDoTest = false;

    private String iDestination = null;

    public TestFTPClient2() {
        this("Test scenario for the FTPClient class.");
    }

    public TestFTPClient2(String aName) {
        super(aName);

        Properties p = super.getPropertiesFile("FTPClient.properties");
        int liDoTest = Integer.parseInt(p.getProperty("performTest").trim());
        this.iDoTest = (liDoTest>0?true:false);

        if(iDoTest) {
            this.iServer = p.getProperty("server").trim();
            this.iUser = p.getProperty("user").trim();
            this.iPassword = p.getProperty("password").trim();
            this.iDestination = p.getProperty("destination").trim();
            if(!iDestination.endsWith("/")) {
                iDestination = iDestination + "/";
            }
        }
    }

    /**
     * This method test the sending of a single text file.
     */
    public void testSendingTextFile() {
        if(iDoTest) {
            try {
                // First get a textfile to send.
                String file = super.getFullFilePath("FTPClient.properties");

                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, this.iPassword);
                ftpClient.sendTextFile(file);

                // Now find the file on the server.
                File f = new File(iDestination + new File(file).getName());
                // See if it exists.
                Assert.assertTrue(f.exists());
                // Compare it line by line.
                BufferedReader brOriginal = new BufferedReader(new FileReader(file));
                BufferedReader brControl = new BufferedReader(new FileReader(f));

                String line = null;
                while((line = brOriginal.readLine()) != null) {
                    Assert.assertEquals(line, brControl.readLine());
                }
                // Now we should also get a 'null' from the destination file.
                Assert.assertTrue(brControl.readLine() == null);

                // Close readers.
                brControl.close();
                brOriginal.close();

                // Delete the destination file.
                while(!f.delete()) {
                }
            } catch(IOException ioe) {
                logger.error(ioe.getMessage(), ioe);
                fail("IOException thrown when attempting to send a single text file: " + ioe.getMessage());
            }
        }
    }

    /**
     * This method test the sending of a single binary file.
     */
    public void testSendingBinaryFile() {
        if(iDoTest) {
            try {
                // First get a binary file to send.
                String file = super.getFullFilePath("testFile.jpg");

                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, this.iPassword);
                ftpClient.sendBinaryFile(file);

                // Now find the file on the server.
                File f = new File(iDestination + new File(file).getName());
                // See if it exists.
                Assert.assertTrue(f.exists());
                // Compare it byte by byte.
                BufferedInputStream biOriginal = new BufferedInputStream(new FileInputStream(file));
                BufferedInputStream biControl = new BufferedInputStream(new FileInputStream(f));

                int current = -1;
                while((current = biOriginal.read()) != -1) {
                    Assert.assertEquals(current, biControl.read());
                }
                // Now we should also get an 'EOF' byte from the destination file.
                Assert.assertTrue(biControl.read() == -1);

                // Close readers.
                biControl.close();
                biOriginal.close();

                // Delete the destination file.
                while(!f.delete()) {
                }
            } catch(IOException ioe) {
                logger.error(ioe.getMessage(), ioe);
                fail("IOException thrown when attempting to send a single binary file: " + ioe.getMessage());
            }
        }
    }

    /**
     * This method test the sending of a few text files.
     */
    public void testSendingTextFiles() {
        if(iDoTest) {
            try {
                // Get two files to send.
                String[] files = new String[2];
                files[0] = super.getFullFilePath("FTPClient.properties");
                files[1] = super.getFullFilePath("enzymes.txt");

                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, this.iPassword);
                ftpClient.sendFiles(files, false);

                for(int i=0;i<files.length;i++) {
                    // Now find the file on the server.
                    File f = new File(iDestination + new File(files[i]).getName());
                    // See if it exists.
                    Assert.assertTrue(f.exists());
                    // Compare it line by line.
                    BufferedReader brOriginal = new BufferedReader(new FileReader(files[i]));
                    BufferedReader brControl = new BufferedReader(new FileReader(f));

                    String line = null;
                    while((line = brOriginal.readLine()) != null) {
                        Assert.assertEquals(line, brControl.readLine());
                    }
                    // Now we should also get a 'null' from the destination file.
                    Assert.assertTrue(brControl.readLine() == null);

                    // Close readers.
                    brControl.close();
                    brOriginal.close();

                    // Delete the destination file.
                    while(!f.delete()) {
                    }
                }
            } catch(IOException ioe) {
                logger.error(ioe.getMessage(), ioe);
                fail("IOException thrown when attempting to send a few text files: " + ioe.getMessage());
            }
        }
    }

    /**
     * This method test the sending of a few binary files.
     */
    public void testSendingBinaryFiles() {
        if(iDoTest) {
            try {
                // Get two files to send.
                String[] files = new String[2];
                files[0] = super.getFullFilePath("testFile.jpg");
                files[1] = super.getFullFilePath("TestMonitor.zip");

                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, this.iPassword);
                ftpClient.sendFiles(files, true);

                for(int i=0;i<files.length;i++) {
                    // Now find the file on the server.
                    File f = new File(iDestination + new File(files[i]).getName());
                    // See if it exists.
                    Assert.assertTrue(f.exists());
                    // Compare it byte by byte.
                    BufferedInputStream biOriginal = new BufferedInputStream(new FileInputStream(files[i]));
                    BufferedInputStream biControl = new BufferedInputStream(new FileInputStream(f));

                    int current = -1;
                    while((current = biOriginal.read()) != -1) {
                        Assert.assertEquals(current, biControl.read());
                    }
                    // Now we should also get an 'EOF' byte from the destination file.
                    Assert.assertTrue(biControl.read() == -1);

                    // Close readers.
                    biControl.close();
                    biOriginal.close();

                    // Delete the destination file.
                    while(!f.delete()) {
                    }
                }
            } catch(IOException ioe) {
                logger.error(ioe.getMessage(), ioe);
                fail("IOException thrown when attempting to send a few binary files: " + ioe.getMessage());
            }
        }
    }

    /**
     * This method test the sending of a group of mixed files.
     */
    public void testSendingMixedFiles() {
        if(iDoTest) {
            try {
                // We get a textfile and a binary file.
                String[] files = new String[2];
                files[0] = super.getFullFilePath("testFile.jpg");
                files[1] = super.getFullFilePath("FTPClient.properties");

                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, this.iPassword);
                ftpClient.sendFiles(files, new boolean[]{true, false});

                // Check the binary file.
                // Now find the file on the server.
                File f = new File(iDestination + new File(files[0]).getName());
                // See if it exists.
                Assert.assertTrue(f.exists());
                // Compare it byte by byte.
                BufferedInputStream biOriginal = new BufferedInputStream(new FileInputStream(files[0]));
                BufferedInputStream biControl = new BufferedInputStream(new FileInputStream(f));

                int current = -1;
                while((current = biOriginal.read()) != -1) {
                    Assert.assertEquals(current, biControl.read());
                }
                // Now we should also get an 'EOF' byte from the destination file.
                Assert.assertTrue(biControl.read() == -1);

                // Close readers.
                biControl.close();
                biOriginal.close();

                // Delete the destination file.
                while(!f.delete()) {
                }

                // Check the text file.
                // Find the file on the server.
                f = new File(iDestination + new File(files[1]).getName());
                // See if it exists.
                Assert.assertTrue(f.exists());
                // Compare it line by line.
                BufferedReader brOriginal = new BufferedReader(new FileReader(files[1]));
                BufferedReader brControl = new BufferedReader(new FileReader(f));

                String line = null;
                while((line = brOriginal.readLine()) != null) {
                    Assert.assertEquals(line, brControl.readLine());
                }
                // Now we should also get a 'null' from the destination file.
                Assert.assertTrue(brControl.readLine() == null);

                // Close readers.
                brControl.close();
                brOriginal.close();

                // Delete the destination file.
                while(!f.delete()) {
                }
            } catch(IOException ioe) {
                logger.error(ioe.getMessage(), ioe);
                fail("IOException thrown when attempting to send mixed files: " + ioe.getMessage());
            }
        }
    }

    /**
     * This method test the test connection method.
     */
    public void testTest() {
        if(iDoTest) {
            // This one should work!
            try {
                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, this.iPassword);
                ftpClient.testFTPConnection();
            } catch(IOException ioe) {
                fail("IOException thrown when testing the connection with correct data: " + ioe.getMessage());
            }

            // These should fail.
            try {
                FTPClient2 ftpClient = new FTPClient2("I_DO_NOT_EXIST_" + this.iServer, this.iUser, this.iPassword);
                ftpClient.testFTPConnection();
                fail("IOException NOT thrown when testing the connection to a non-existing server!");
            } catch(IOException ioe) {
                // We want this to happen.
            }

            try {
                FTPClient2 ftpClient = new FTPClient2(this.iServer, "I_DO_NOT_EXIST_" + this.iUser, this.iPassword);
                ftpClient.testFTPConnection();
                fail("IOException NOT thrown when testing the connection to a non-existing user!");
            } catch(IOException ioe) {
                // We want this to happen.
            }

            try {
                FTPClient2 ftpClient = new FTPClient2(this.iServer, this.iUser, "I_DO_NOT_EXIST_" + this.iPassword);
                ftpClient.testFTPConnection();
                fail("IOException NOT thrown when testing the connection to a non-existing password!");
            } catch(IOException ioe) {
                // We want this to happen.
            }
        }
    }
}
