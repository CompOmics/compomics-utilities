package com.compomics.util.experiment.biology.ions;

/**
 * This class contains convenience methods for the handling of charges.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class Charge {

    /**
     * Empty default constructor.
     */
    public Charge() {
    }

    /**
     * Returns a string representing the charge. For example 2+.
     *
     * @param value the value of the charge
     *
     * @return charge as a string
     */
    public static String toString(int value) {

        if (value == 0) {
            return "0";
        }

        return value > 0 ? value + "+" : value + "-";
    }
}
