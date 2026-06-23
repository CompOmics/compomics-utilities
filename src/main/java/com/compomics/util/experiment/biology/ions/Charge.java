package com.compomics.util.experiment.biology.ions;

/**
 * This class contains convenience methods for the handling of charges.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class Charge {

    /**
     * Legacy charge value field.
     */
    public int value;

    /**
     * Empty default constructor.
     */
    public Charge() {
    }

    /**
     * Constructor with a charge value.
     *
     * @param value the charge value
     */
    public Charge(int value) {
        this.value = value;
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

    /**
     *
     * Returns the charge as a string of + or -. One for each charge. A charge
     * of +1 however returns the empty string.
     *
     * @param value the value of the charge
     *
     * @return the charge as a string of +
     */
    public static String getChargeAsFormattedString(int value) {

        if (value == 1) {
            return "";
        }

        int absValue = Math.abs(value);

        StringBuilder chargeAsString = new StringBuilder(absValue);

        for (int i = 0; i < absValue; i++) {

            if (value > 0) {

                chargeAsString.append('+');

            } else {

                chargeAsString.append('-');

            }
        }

        return chargeAsString.toString();

    }

}
