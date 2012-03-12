package com.compomics.util.pride.prideobjects;

import com.compomics.util.pride.CvTerm;
import com.compomics.util.pride.PrideObject;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * An object for storing Instrument details.
 *
 * @author Harald Barsnes
 */
public class Instrument implements PrideObject, Serializable {

    /**
     * serialization number for backward compatibility
     */
    static final long serialVersionUID = -8802861658166703745L;
    /**
     * The instrument name.
     */
    private String name;
    /**
     * The list of analyzer CV terms.
     */
    private ArrayList<CvTerm> cvTerms;
    /**
     * The instrument source.
     */
    private CvTerm source;
    /**
     * The instrument detector.
     */
    private CvTerm detector;

    /**
     * Create a new Instrument object.
     *
     * @param name
     * @param source
     * @param detector
     * @param cvTerms
     */
    public Instrument(String name, CvTerm source, CvTerm detector, ArrayList<CvTerm> cvTerms) {
        this.name = name;
        this.source = source;
        this.detector = detector;
        this.cvTerms = cvTerms;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the cvTerms
     */
    public ArrayList<CvTerm> getCvTerms() {
        return cvTerms;
    }

    /**
     * @param cvTerms the cvTerms to set
     */
    public void setCvTerms(ArrayList<CvTerm> cvTerms) {
        this.cvTerms = cvTerms;
    }

    /**
     * @return the source
     */
    public CvTerm getSource() {
        return source;
    }

    /**
     * @param source the source to set
     */
    public void setSource(CvTerm source) {
        this.source = source;
    }

    /**
     * @return the detector
     */
    public CvTerm getDetector() {
        return detector;
    }

    /**
     * @param detector the detector to set
     */
    public void setDetector(CvTerm detector) {
        this.detector = detector;
    }

    /**
     * Returns a list of predefined instruments
     *
     * @return a list of predefined instruments
     */
    public static ArrayList<Instrument> getDefaultInstruments() {
        ArrayList<Instrument> result = new ArrayList<Instrument>();
        result.add(new Instrument("Bruker Ultraflex",
                new CvTerm("MS", "MS:1000075", "Matrix-assisted Laser Desorption Ionization", null),
                new CvTerm("MS", "MS:1000111", "Electron Multiplier Tube", null),
                new ArrayList<CvTerm>(Arrays.asList(
                new CvTerm("MS", "MS:1000202", "Bruker Daltonics ultraFlex TOF/TOF MS", null)))));
        result.add(new Instrument("LCQ Duo",
                new CvTerm("MS", "MS:1000073", "Electrospray Ionization", null),
                new CvTerm("MS", "MS:1000111", "Electron Multiplier Tube", null),
                new ArrayList<CvTerm>(Arrays.asList(
                new CvTerm("MS", "MS:1000264", "Ion Trap", null)))));
        result.add(new Instrument("LTQ-Orbitrap",
                new CvTerm("MS", "MS:1000073", "Electrospray Ionization", null),
                new CvTerm("MS", "MS:1000111", "Electron Multiplier Tube", null),
                new ArrayList<CvTerm>(Arrays.asList(
                new CvTerm("MS", "MS:1000449", "LTQ Orbitrap", null)))));
        result.add(new Instrument("QSTAR Pulsar I",
                new CvTerm("MS", "MS:1000075", "Matrix-assisted Laser Desorption Ionization", null),
                new CvTerm("MS", "MS:1000108", "Conversion Dynode Electron Multiplier", null),
                new ArrayList<CvTerm>(Arrays.asList(
                new CvTerm("MS", "MS:1000084", "Time-of-flight", null)))));
        result.add(new Instrument("QToF Global",
                new CvTerm("MS", "MS:1000114", "Electrospray Ionization", null),
                new CvTerm("MS", "MS:1000111", "Microchannel Plate Detector", null),
                new ArrayList<CvTerm>(Arrays.asList(
                new CvTerm("MS", "MS:1000081", "Quadrupole", null),
                new CvTerm("MS", "MS:1000084", "Time-of-flight", null)))));
        result.add(new Instrument("QToF-Ultima",
                new CvTerm("MS", "MS:1000073", "Electrospray Ionization", null),
                new CvTerm("MS", "MS:1000111", "Electron Multiplier Tube", null),
                new ArrayList<CvTerm>(Arrays.asList(
                new CvTerm("MS", "MS:1000188", "Q-Tof micro", null)))));
        return result;
    }

    public String getFileName() {
        return name;
    }
}
