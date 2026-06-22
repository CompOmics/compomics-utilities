package com.compomics.util.exceptions;

import java.util.HashSet;

/**
 * Interface for the general handling of exceptions.
 *
 * @author Marc Vaudel
 */
public abstract class ExceptionHandler {

    /**
     * Empty default constructor
     */
    public ExceptionHandler() {
    }

    /**
     * List of caught exceptions.
     */
    protected HashSet<String> exceptionCaught = new HashSet<>();
    /**
     * If true exceptions will be ignored.
     */
    protected boolean ignoreExceptions = false;

    /**
     * Catches an exception and informs the user.
     *
     * @param e the exception caught
     */
    public synchronized void catchException(Exception e) {
        
        if (!ignoreExceptions && !exceptionCaught.contains(getExceptionType(e))) {

            e.printStackTrace();
            exceptionCaught.add(getExceptionType(e));

            // @TODO: remove once the underlying Nimbus look and feel bug is fixed.
            // On recent JDKs the Nimbus look and feel can throw a benign
            // ClassCastException ("ColorUIResource cannot be cast to Boolean" in
            // NimbusStyle.isOpaque) while building chart popup menus. It does not
            // affect functionality, so it is logged above but not shown to the user.
            if (isBenignLookAndFeelException(e)) {
                return;
            }

            notifyUser(e);

        }
    }

    /**
     * Indicates whether the given exception is the known benign look and feel
     * ClassCastException thrown while rendering (e.g. "ColorUIResource cannot be
     * cast to Boolean" in NimbusStyle). Such exceptions do not affect
     * functionality and should not be reported to the user.
     *
     * @param e the exception to inspect
     *
     * @return true if the exception is a benign look and feel rendering exception
     */
    private static boolean isBenignLookAndFeelException(Exception e) {

        if (!(e instanceof ClassCastException)) {
            return false;
        }

        for (StackTraceElement element : e.getStackTrace()) {

            String className = element.getClassName();

            if (className.startsWith("javax.swing.plaf.nimbus.")
                    || className.startsWith("javax.swing.plaf.synth.")) {
                return true;
            }
        }

        return false;
    }

    /**
     * Notifies the user that an exception was caught.
     *
     * @param e the exception to catch
     */
    protected abstract void notifyUser(Exception e);

    /**
     * Returns the exception type.
     *
     * @param e the exception to get the type fro
     * @return the exception type as a string
     */
    public static String getExceptionType(Exception e) {
        if (e.getLocalizedMessage() == null) {
            return "null pointer";
        } else if (e.getLocalizedMessage().startsWith("Protein not found")) {
            return "Protein not found";
        } else if (e.getLocalizedMessage().startsWith("Error while loading")
                || e.getLocalizedMessage().startsWith("Error while writing")) {
            return "Serialization";
        } else if (e.getLocalizedMessage().startsWith("Two modifications found")) {
            return "Two modifications on same site";
        } else {
            return e.getLocalizedMessage();
        }
    }

    /**
     * Sets whether exceptions should be ignored.
     *
     * @param ignoreExceptions if true exceptions will be ignored
     */
    public void setIgnoreExceptions(boolean ignoreExceptions) {
        this.ignoreExceptions = ignoreExceptions;
    }
}
