package com.compomics.util.experiment.identification.matches;

import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.ions.impl.TagFragmentIon;
import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;
import com.compomics.util.experiment.biology.atoms.Atom;
import com.compomics.util.experiment.biology.ions.Charge;
import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.identification.spectrum_annotation.IonMatchKeysCache;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.pride.CvTerm;

/**
 * This class represents the assignment of a peak to a theoretical ion.
 *
 * @author Marc Vaudel
 */
public class IonMatch extends ExperimentObject {

    /**
     * Empty default constructor
     */
    public IonMatch() {
    }
    
    /**
     * The matched peak m/z.
     */
    public double peakMz;
    /**
     * The matched peak intensity.
     */
    public double peakIntensity;
    /**
     * The matching ion.
     */
    public Ion ion;
    /**
     * The inferred charge of the ion.
     */
    public int charge;

    /**
     * Constructor for an ion match.
     *
     * @param peakMz The matched peak m/z.
     * @param peakIntensity The matched peak intensity.
     * @param ion The theoretic ion.
     * @param charge The inferred charge of the ion.
     */
    public IonMatch(
            double peakMz,
            double peakIntensity,
            Ion ion,
            int charge
    ) {
        this.peakMz = peakMz;
        this.peakIntensity = peakIntensity;
        this.ion = ion;
        this.charge = charge;
    }

    /**
     * Get the absolute matching error in Da.
     *
     * @return the absolute matching error
     */
    public double getAbsoluteError() {
        readDBMode();
        double theoreticMz = ion.getTheoreticMz(charge);
        return peakMz - theoreticMz;
    }

    /**
     * Get the absolute matching error in Da after isotope removal.
     *
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     *
     * @return the absolute matching error
     */
    public double getAbsoluteError(
            int minIsotope,
            int maxIsotope
    ) {
        readDBMode();
        double theoreticMz = ion.getTheoreticMz(charge);
        double measuredMz = peakMz;
        measuredMz -= getIsotopeNumber(minIsotope, maxIsotope) * Atom.C.getDifferenceToMonoisotopic(1) / charge;
        return measuredMz - theoreticMz;
    }

    /**
     * Get the relative m/z matching error in ppm.
     *
     * @return the relative matching error
     */
    public double getRelativeError() {
        readDBMode();
        double theoreticMz = ion.getTheoreticMz(charge);
        return ((peakMz - theoreticMz) * 1000000) / theoreticMz;
    }

    /**
     * Get the relative m/z matching error in ppm after isotope removal.
     *
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     *
     * @return the relative matching error
     */
    public double getRelativeError(
            int minIsotope,
            int maxIsotope
    ) {
        readDBMode();
        double theoreticMz = ion.getTheoreticMz(charge);
        double measuredMz = peakMz;
        measuredMz -= getIsotopeNumber(minIsotope, maxIsotope) * Atom.C.getDifferenceToMonoisotopic(1) / charge;
        return ((measuredMz - theoreticMz) * 1000000) / theoreticMz;
    }

    /**
     * Returns the distance in number of neutrons between the experimental mass
     * and theoretic mass, image of the isotope number: 1 typically indicates
     * C13 isotope.
     *
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     *
     * @return the distance in number of neutrons between the experimental mass
     * and theoretic mass
     */
    public int getIsotopeNumber(
            int minIsotope,
            int maxIsotope
    ) {
        readDBMode();
        double experimentalMass = peakMz * charge - charge * ElementaryIon.proton.getTheoreticMass();
        double result = (experimentalMass - ion.getTheoreticMass()) / Atom.C.getDifferenceToMonoisotopic(1);
        return Math.min(Math.max((int) Math.round(result), minIsotope), maxIsotope);
    }

    /**
     * Returns the error.
     *
     * @param isPpm a boolean indicating whether the error should be retrieved
     * in ppm (true) or in Dalton (false)
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     *
     * @return the match m/z error
     */
    public double getError(
            boolean isPpm,
            int minIsotope,
            int maxIsotope
    ) {

        return isPpm ? getRelativeError(minIsotope, maxIsotope) : getAbsoluteError(minIsotope, maxIsotope);

    }

    /**
     * Returns the error.
     *
     * @param isPpm a boolean indicating whether the error should be retrieved
     * in ppm (true) or in Dalton (false)
     *
     * @return the match m/z error
     */
    public double getError(
            boolean isPpm
    ) {
        readDBMode();

        return isPpm ? getRelativeError() : getAbsoluteError();

    }

    /**
     * Returns the annotation to use for the ion match as a String.
     *
     * @return the annotation to use for the given ion match
     */
    public String getPeakAnnotation() {
        readDBMode();
        return getPeakAnnotation(
                false,
                ion,
                charge
        );
    }

    /**
     * Returns the annotation to use for a given ion and charge as a String.
     *
     * @param ion the given ion
     * @param charge the given charge
     *
     * @return the annotation to use for the given ion match
     */
    public static String getPeakAnnotation(
            Ion ion,
            int charge
    ) {
        return getPeakAnnotation(
                false,
                ion,
                charge
        );
    }

    /**
     * Returns the key for the ion match uniquely representing a peak
     * annotation.
     *
     * @param ion the ion matched
     * @param charge the charge
     *
     * @return the key for the ion match
     */
    public static String getMatchKey(
            Ion ion,
            int charge
    ) {
        return getMatchKey(
                ion,
                charge,
                null
        );
    }

    /**
     * Returns the key for the ion match uniquely representing a peak
     * annotation. If a cache is given it will be used to store keys, ignored if
     * null.
     *
     * @param ion the ion matched
     * @param charge the charge
     * @param ionMatchKeysCache a cache for the ion match keys
     *
     * @return the key for the ion match
     */
    public static String getMatchKey(
            Ion ion,
            int charge,
            IonMatchKeysCache ionMatchKeysCache
    ) {

        if (ionMatchKeysCache != null) {

            return ionMatchKeysCache.getMatchKey(ion, charge);

        }

        Ion.IonType ionType = ion.getType();
        int ionTypeIndex = ionType.index;
        int ionSubType = ion.getSubType();
        int fragmentIonNumber;

        if (ionType == Ion.IonType.PEPTIDE_FRAGMENT_ION) {

            PeptideFragmentIon fragmentIon = ((PeptideFragmentIon) ion);
            fragmentIonNumber = fragmentIon.getNumber();

        } else if (ionType == Ion.IonType.TAG_FRAGMENT_ION) {

            TagFragmentIon tagFragmentIon = ((TagFragmentIon) ion);
            fragmentIonNumber = tagFragmentIon.getNumber();

        } else {

            fragmentIonNumber = 0;

        }

        String neutralLossesAsString = ion.getNeutralLossesAsString();
        String key = getMatchKey(
                ionTypeIndex,
                ionSubType,
                fragmentIonNumber,
                neutralLossesAsString,
                charge
        );

        return key;

    }

    /**
     * Returns the key based on the different attributes of a match.
     *
     * @param ionTypeIndex the index of the ion type
     * @param ionSubType the index of the ion subtype
     * @param fragmentIonNumber the number of the ion, 0 if none
     * @param neutralLossesAsString the neutral losses as a string
     * @param charge the charge
     *
     * @return the key for the ion match
     */
    public static String getMatchKey(
            int ionTypeIndex,
            int ionSubType,
            int fragmentIonNumber,
            String neutralLossesAsString,
            int charge
    ) {

        StringBuilder stringBuilder = new StringBuilder(8);

        stringBuilder
                .append(ionTypeIndex)
                .append("_")
                .append(ionSubType)
                .append("_")
                .append(fragmentIonNumber)
                .append("_")
                .append(neutralLossesAsString)
                .append("_")
                .append(charge);

        return stringBuilder.toString();

    }

    /**
     * Returns the annotation to use for a given ion and charge as a String.
     *
     * @param html if true, returns the annotation as HTML with subscripts tags
     * @param ion the given ion
     * @param charge the given charge
     *
     * @return the annotation to use for the given ion match
     */
    public static String getPeakAnnotation(
            boolean html,
            Ion ion,
            int charge
    ) {

        StringBuilder result = new StringBuilder();

        switch (ion.getType()) {
            case PEPTIDE_FRAGMENT_ION:
                if (html) {
                    result.append("<html>");
                }
                result.append(ion.getSubTypeAsString());

                // add fragment ion number
                PeptideFragmentIon fragmentIon = ((PeptideFragmentIon) ion);
                if (html) {
                    result.append("<sub>").append(fragmentIon.getNumber()).append("</sub>");
                } else {
                    result.append(fragmentIon.getNumber());
                }

                // add charge
                result.append(Charge.toString(charge));

                // add any neutral losses
                if (html) {
                    String neutralLoss = ion.getNeutralLossesAsString();

                    for (int i = 0; i < neutralLoss.length(); i++) {
                        if (Character.isDigit(neutralLoss.charAt(i))) {
                            result.append("<sub>").append(neutralLoss.charAt(i)).append("</sub>");
                        } else {
                            result.append(neutralLoss.charAt(i));
                        }
                    }
                } else {
                    result.append(ion.getNeutralLossesAsString());
                }
                if (html) {
                    result.append("</html>");
                }
                return result.toString();

            case TAG_FRAGMENT_ION:
                TagFragmentIon tagFragmentIon = (TagFragmentIon) ion;

                if (html) {
                    result.append("<html>");
                }
                // add type
                result.append(ion.getSubTypeAsString());

                // add fragment ion number
                if (html) {
                    result.append("<sub>").append(tagFragmentIon.getSubNumber()).append("</sub>");
                } else {
                    result.append(tagFragmentIon.getSubNumber());
                }

                // add charge
                result.append(Charge.toString(charge));

                // add any neutral losses
                if (html) {
                    String neutralLoss = ion.getNeutralLossesAsString();

                    for (int i = 0; i < neutralLoss.length(); i++) {
                        if (Character.isDigit(neutralLoss.charAt(i))) {
                            result.append("<sub>").append(neutralLoss.charAt(i)).append("</sub>");
                        } else {
                            result.append(neutralLoss.charAt(i));
                        }
                    }
                } else {
                    result.append(ion.getNeutralLossesAsString());
                }

                if (html) {
                    result.append("</html>");
                }
                return result.toString();

            case PRECURSOR_ION:
                if (html) {
                    result.append("<html>");
                }
                result.append(ion.getSubTypeAsString()).append("-");

                // add charge
                result.append(Charge.toString(charge));

                // add any neutral losses
                String neutralLoss = ion.getNeutralLossesAsString();
                if (html) {
                    for (int i = 0; i < neutralLoss.length(); i++) {
                        if (Character.isDigit(neutralLoss.charAt(i))) {
                            result.append("<sub>").append(neutralLoss.charAt(i)).append("</sub>");
                        } else {
                            result.append(neutralLoss.charAt(i));
                        }
                    }
                } else {
                    result.append(neutralLoss);
                }
                if (html) {
                    result.append("</html>");
                }
                return result.toString();

            default:
                if (html) {
                    result.append("<html>");
                }
                result.append(ion.getName());
                if (html) {
                    result.append("</html>");
                }
                return result.toString();

        }
    }

    /**
     * Returns the annotation to use for the given ion match as a String.
     *
     * @param html if true, returns the annotation as HTML with subscripts tags
     * @return the annotation to use for the given ion match
     */
    public String getPeakAnnotation(
            boolean html
    ) {
        readDBMode();
        return getPeakAnnotation(
                html,
                ion,
                charge
        );
    }

    /**
     * Returns the pride CV term for the ion match m/z.
     *
     * @return the pride CV term for the ion match m/z
     */
    public CvTerm getMZPrideCvTerm() {
        readDBMode();
        return new CvTerm(
                "PRIDE",
                "PRIDE:0000188",
                "product ion m/z",
                Double.toString(peakMz)
        );
    }

    /**
     * Returns the pride CV term for the ion match intensity.
     *
     * @return the pride CV term for the ion match intensity
     */
    public CvTerm getIntensityPrideCvTerm() {
        readDBMode();
        return new CvTerm(
                "PRIDE", 
                "PRIDE:0000189", 
                "product ion intensity", 
                Double.toString(peakIntensity)
        );
    }

    /**
     * Returns the pride CV term for the ion match error.
     *
     * @param minIsotope the minimal isotope
     * @param maxIsotope the maximal isotope
     *
     * @return the pride CV term for the ion match error
     */
    public CvTerm getIonMassErrorPrideCvTerm(
            int minIsotope, 
            int maxIsotope
    ) {
        readDBMode();
        return new CvTerm(
                "PRIDE", 
                "PRIDE:0000190", 
                "product ion mass error", 
                Double.toString(
                        getAbsoluteError(
                                minIsotope, 
                                maxIsotope
                        )
                )
        );
    }

    /**
     * Returns the pride CV term for the ion match charge.
     *
     * @return the pride CV term for the ion match charge
     */
    public CvTerm getChargePrideCvTerm() {
        readDBMode();
        return new CvTerm(
                "PRIDE", 
                "PRIDE:0000204", 
                "product ion charge", 
                Integer.toString(charge)
        );
    }

    /**
     * Enum of the supported error types.
     */
    public enum MzErrorType {

        Absolute(
                "Absolute", 
                "Absolute error", 
                "m/z"
        ),
        RelativePpm(
                "Relative (ppm)", 
                "Relative error in ppm", 
                "ppm"
        ),
        Statistical(
                "Statistical", 
                "Probability to reach this error according to the error distribution", 
                "%p"
        );
        /**
         * The name of the error type.
         */
        public final String name;
        /**
         * The description of the error type.
         */
        public final String description;
        /**
         * The unit to use
         */
        public final String unit;

        /**
         * Constructor.
         *
         * @param name the name of the error type
         * @param description the description of the error type
         * @param unit the unit to use
         */
        private MzErrorType(
                String name, 
                String description, 
                String unit
        ) {
            this.name = name;
            this.description = description;
            this.unit = unit;
        }

        /**
         * Returns the error type corresponding to the given index. Error types
         * are indexed according to the values() method. Null if not found.
         *
         * @param index the index of the error type in the values() method
         *
         * @return the corresponding error type
         */
        public static MzErrorType getMzErrorType(
                int index
        ) {
        
            MzErrorType[] values = MzErrorType.values();
            
            if (index >= 0 && index < values.length) {
            
                return values[index];
            
            }
            
            return null;
        
        }
    }
}
