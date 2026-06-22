package com.compomics.util.experiment.mass_spectrometry.spectra;

import com.compomics.util.experiment.biology.ions.Charge;
import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;
import com.compomics.util.experiment.personalization.ExperimentObject;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;

/**
 * This class models a precursor.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class Precursor extends ExperimentObject {

    /**
     * The retention time when the precursor was isolated.
     */
    public final double rt;
    /**
     * The measured m/z of the precursor.
     */
    public final double mz;
    /**
     * The measured intensity of the precursor.
     */
    public final double intensity;
    /**
     * The charge of the precursor.
     */
    public final int[] possibleCharges;

    /**
     * Empty default constructor.
     */
    public Precursor() {

        this.rt = Double.NaN;
        this.mz = Double.NaN;
        this.intensity = Double.NaN;
        this.possibleCharges = new int[0];

    }

    /**
     * Constructor for the precursor.
     *
     * @param rt the retention time
     * @param mz the m/z
     * @param possibleCharges the possible charges
     */
    public Precursor(
            double rt,
            double mz,
            int[] possibleCharges
    ) {

        this.rt = rt;
        this.mz = mz;
        this.intensity = Double.NaN;
        this.possibleCharges = new int[possibleCharges.length];

        System.arraycopy(possibleCharges, 0, this.possibleCharges, 0, possibleCharges.length);

    }

    /**
     * Constructor with retention time window.
     *
     * @param rt the retention time
     * @param mz the m/z
     * @param intensity the intensity
     * @param possibleCharges the possible charges
     */
    public Precursor(
            double rt,
            double mz,
            double intensity,
            int[] possibleCharges
    ) {

        this.rt = rt;
        this.mz = mz;
        this.intensity = intensity;
        this.possibleCharges = new int[possibleCharges.length];
        System.arraycopy(possibleCharges, 0, this.possibleCharges, 0, possibleCharges.length);

    }

    /**
     * Constructor with retention time window and no reference retention time.
     *
     * @param mz the m/z
     * @param intensity the intensity
     * @param possibleCharges the possible charges
     * @param rtMin the minimum of the RT window in seconds
     * @param rtMax the maximum of the RT window in seconds
     */
    public Precursor(
            double mz,
            double intensity,
            int[] possibleCharges,
            double rtMin,
            double rtMax
    ) {

        this.rt = (rtMin + rtMax) / 2;
        this.mz = mz;
        this.intensity = intensity;
        this.possibleCharges = new int[possibleCharges.length];

        System.arraycopy(possibleCharges, 0, this.possibleCharges, 0, possibleCharges.length);

    }

    /**
     * Returns the retention time in minutes.
     *
     * @return the retention time in minutes
     */
    public double getRtInMinutes() {

        return rt / 60;

    }

    /**
     * Returns the retention time in seconds.
     *
     * @return the retention time in seconds
     */
    public double getRt() {
        return rt;
    }

    /**
     * Returns the measured precursor m/z.
     *
     * @return the measured precursor m/z
     */
    public double getMz() {
        return mz;
    }

    public double getIntensity() {
        return intensity;
    }

    public ArrayList<Charge> getPossibleCharges() {
        ArrayList<Charge> result = new ArrayList<>(possibleCharges.length);
        for (int possibleCharge : possibleCharges) {
            result.add(new Charge(possibleCharge));
        }
        return result;
    }

    /**
     * Returns the possible charges as a string.
     *
     * @return the possible charges as a string
     */
    public String getPossibleChargesAsString() {

        return Arrays.stream(possibleCharges)
                .mapToObj(
                        charge -> Charge.toString(charge)
                )
                .collect(Collectors.joining(", "));

    }

    /**
     * Returns the mass of the precursor with the given charge.
     *
     * @param chargeValue the value of the charge
     *
     * @return the mass of the precursor with the given charge
     */
    public double getMass(
            int chargeValue
    ) {

        return mz * chargeValue - chargeValue * ElementaryIon.proton.getTheoreticMass();

    }

    /**
     * Returns a boolean indicating whether the precursor is identical to the
     * other precursor. m/z, rt, and intensities values must have exact same
     * double values. Charges must be identical and in the same order.
     *
     * @param otherPrecursor The other precursor.
     *
     * @return A boolean indicating whether the precursor is identical to the
     * other precursor.
     */
    public boolean isSameAs(
            Precursor otherPrecursor
    ) {

        if (possibleCharges.length != otherPrecursor.possibleCharges.length) {
            return false;
        }

        for (int i = 0; i < possibleCharges.length; i++) {

            if (possibleCharges[i] != otherPrecursor.possibleCharges[i]) {

                return false;

            }
        }

        if (Double.isNaN(mz) && !Double.isNaN(otherPrecursor.mz)
                || !Double.isNaN(mz) && Double.isNaN(otherPrecursor.mz)
                || !Double.isNaN(mz) && !Double.isNaN(otherPrecursor.mz) && mz != otherPrecursor.mz) {

            return false;

        }

        if (Double.isNaN(intensity) && !Double.isNaN(otherPrecursor.intensity)
                || !Double.isNaN(intensity) && Double.isNaN(otherPrecursor.intensity)
                || !Double.isNaN(intensity) && !Double.isNaN(otherPrecursor.intensity) && intensity != otherPrecursor.intensity) {

            return false;

        }

        if (Double.isNaN(rt) && !Double.isNaN(otherPrecursor.rt)
                || !Double.isNaN(rt) && Double.isNaN(otherPrecursor.rt)
                || !Double.isNaN(rt) && !Double.isNaN(otherPrecursor.rt) && rt != otherPrecursor.rt) {

            return false;

        }

        return true;

    }
}
