package com.compomics.util.experiment.identification.ptm;

import java.util.ArrayList;

/**
 * An enum of the PTM scores.
 *
 * @author Marc Vaudel
 */
public enum PtmScore {

    /**
     * The A score.
     */
    AScore(0, "A-score"),
    /**
     * The PhosphoRS score.
     */
    PhosphoRS(1, "PhosphoRS");
    /**
     * Score id number.
     */
    private int id;
    /**
     * Score name.
     */
    private String name;

    /**
     * Constructor.
     *
     * @param id the id number
     * @param name the name
     */
    private PtmScore(int id, String name) {
        this.id = id;
        this.name = name;
    }

    /**
     * Returns the id number of the score.
     *
     * @return the id number of the score
     */
    public int getId() {
        return id;
    }

    /**
     * Returns the name of the score.
     *
     * @return the name of the score
     */
    public String getName() {
        return name;
    }

    /**
     * Returns a list of the implemented scores.
     *
     * @return a list of the implemented scores
     */
    public static ArrayList<PtmScore> getImplementedPtmScores() {
        ArrayList<PtmScore> result = new ArrayList<PtmScore>();
        result.add(AScore);
        result.add(PhosphoRS);
        return result;
    }

    /**
     * Returns the PTM score indexed by the given id.
     *
     * @param id the id number of the PTM score
     * @return the desired PTM score
     */
    public static PtmScore getScore(int id) {
        for (PtmScore ptmScore : getImplementedPtmScores()) {
            if (ptmScore.getId() == id) {
                return ptmScore;
            }
        }
        throw new IllegalArgumentException("Ptm score of id " + id + " not recognized.");
    }

    /**
     * Returns the PTM score of the given name.
     *
     * @param name the name of the score
     * @return the desired PTM score
     */
    public static PtmScore getScore(String name) {
        for (PtmScore ptmScore : getImplementedPtmScores()) {
            if (ptmScore.getName().equals(name)) {
                return ptmScore;
            }
        }
        throw new IllegalArgumentException("Ptm score of name " + name + " not recognized.");
    }

    /**
     * Returns the different implemented scores as list of command line option
     *
     * @return the different implemented scores as list of command line option
     */
    public static String getCommandLineOptions() {
        String result = "";
        for (PtmScore ptmScore : getImplementedPtmScores()) {
            if (!result.equals("")) {
                result += ", ";
            }
            result += ptmScore.getId() + ": " + ptmScore.getName();
        }
        return result;
    }

    /**
     * Returns a list containing the names of the implemented scores.
     *
     * @return a list containing the names of the implemented scores
     */
    public static String[] getScoreNames() {
        ArrayList<PtmScore> scores = getImplementedPtmScores();
        String[] names = new String[scores.size()];
        for (int i = 0; i < scores.size(); i++) {
            names[i] = scores.get(i).getName();
        }
        return names;
    }
}
