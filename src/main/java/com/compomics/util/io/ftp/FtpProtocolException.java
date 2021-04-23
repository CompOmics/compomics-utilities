/*
 * Copyright (C) Lennart Martens
 * 
 * Contact: lennart.martens AT UGent.be (' AT ' to be replaced with '@')
 */

/*
 * Created by IntelliJ IDEA.
 * User: Lennart
 * Date: 29-nov-02
 * Time: 13:54:53
 */
package com.compomics.util.io.ftp;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import java.io.IOException;

/*
 * CVS information:
 *
 * $Revision: 1.3 $
 * $Date: 2007/07/06 09:41:53 $
 */

/**
 * This class handles an FTP protocol exception.
 *
 * @author Lennart
 */
public class FtpProtocolException extends IOException {

    /**
     * Empty default constructor
     */
    public FtpProtocolException() {
    }

    // Class specific log4j logger for FtpProtocolException instances.
    Logger logger = LogManager.getLogger(FtpProtocolException.class);

    /**
     * Creates a new FtpProtocolException from the given string.
     *
     * @param s the FTP protocol exception as a string
     */
    FtpProtocolException(String s) {
        super(s);
    }
}
