package com.compomics.util.parameters.identification.tool_specific;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.gui.parameters.identification.IdentificationAlgorithmParameter;

/**
 * InstaNovo specific parameters.
 *
 * @author CompOmics
 */
public class InstaNovoParameters extends ExperimentObject implements IdentificationAlgorithmParameter {

    /**
     * Version number for deserialization.
     */
    static final long serialVersionUID = -2295564912139753378L;
    /**
     * Default InstaNovo model identifier for v1.2.2.
     */
    public static final String DEFAULT_INSTANOVO_MODEL = "instanovo-v1.2.0";
    /**
     * Default InstaNovo+ model identifier for v1.2.2 refinement.
     */
    public static final String DEFAULT_INSTANOVO_PLUS_MODEL = "instanovoplus-v1.1.0";
    /**
     * Default prediction batch size for desktop SearchGUI runs.
     */
    public static final int DEFAULT_BATCH_SIZE = 16;
    /**
     * The selected InstaNovo model id or path.
     */
    private String instaNovoModel = DEFAULT_INSTANOVO_MODEL;
    /**
     * The selected InstaNovo+ model id or path used for refinement.
     */
    private String instaNovoPlusModel = DEFAULT_INSTANOVO_PLUS_MODEL;
    /**
     * The optional inference configuration path.
     */
    private String configFile = null;
    /**
     * The number of beams.
     */
    private int numberOfBeams = 5;
    /**
     * The prediction batch size. A value below one lets InstaNovo use its
     * configuration default.
     */
    private int batchSize = DEFAULT_BATCH_SIZE;
    /**
     * Whether to use knapsack beam search.
     */
    private boolean useKnapsack = false;
    /**
     * Whether to save all beam search predictions.
     */
    private boolean saveAllPredictions = true;
    /**
     * Whether to force CPU execution.
     */
    private boolean forceCpu = false;

    @Override
    public Advocate getAlgorithm() {
        return Advocate.instanovo;
    }

    @Override
    public boolean equals(IdentificationAlgorithmParameter identificationAlgorithmParameter) {

        if (identificationAlgorithmParameter instanceof InstaNovoParameters) {

            InstaNovoParameters other = (InstaNovoParameters) identificationAlgorithmParameter;

            return safeEquals(instaNovoModel, other.getInstaNovoModel())
                    && safeEquals(instaNovoPlusModel, other.getInstaNovoPlusModel())
                    && safeEquals(configFile, other.getConfigFile())
                    && numberOfBeams == other.getNumberOfBeams()
                    && batchSize == other.getBatchSize()
                    && useKnapsack == other.isUseKnapsack()
                    && saveAllPredictions == other.isSaveAllPredictions()
                    && forceCpu == other.isForceCpu();
        }

        return false;
    }

    @Override
    public String toString(boolean html) {

        String newLine = html ? "<br>" : System.getProperty("line.separator");
        StringBuilder output = new StringBuilder();
        Advocate advocate = getAlgorithm();
        output.append("# ------------------------------------------------------------------");
        output.append(newLine);
        output.append("# ").append(advocate.getName()).append(" Specific Parameters");
        output.append(newLine);
        output.append("# ------------------------------------------------------------------");
        output.append(newLine);
        output.append(newLine);
        output.append("INSTANOVO_MODEL=").append(instaNovoModel).append(newLine);
        output.append("INSTANOVO_PLUS_MODEL=").append(instaNovoPlusModel).append(newLine);
        output.append("CONFIG_FILE=").append(configFile == null ? "" : configFile).append(newLine);
        output.append("NUMBER_OF_BEAMS=").append(numberOfBeams).append(newLine);
        output.append("BATCH_SIZE=").append(batchSize).append(newLine);
        output.append("USE_KNAPSACK=").append(useKnapsack).append(newLine);
        output.append("SAVE_ALL_PREDICTIONS=").append(saveAllPredictions).append(newLine);
        output.append("FORCE_CPU=").append(forceCpu).append(newLine);

        return output.toString();
    }

    /**
     * Returns the selected InstaNovo model.
     *
     * @return the selected InstaNovo model
     */
    public String getInstaNovoModel() {
        return instaNovoModel;
    }

    /**
     * Sets the selected InstaNovo model.
     *
     * @param instaNovoModel the selected InstaNovo model
     */
    public void setInstaNovoModel(String instaNovoModel) {
        this.instaNovoModel = instaNovoModel;
    }

    /**
     * Returns the selected InstaNovo+ model.
     *
     * @return the selected InstaNovo+ model
     */
    public String getInstaNovoPlusModel() {
        return instaNovoPlusModel;
    }

    /**
     * Sets the selected InstaNovo+ model.
     *
     * @param instaNovoPlusModel the selected InstaNovo+ model
     */
    public void setInstaNovoPlusModel(String instaNovoPlusModel) {
        this.instaNovoPlusModel = instaNovoPlusModel;
    }

    /**
     * Returns the optional configuration file.
     *
     * @return the optional configuration file
     */
    public String getConfigFile() {
        return configFile;
    }

    /**
     * Sets the optional configuration file.
     *
     * @param configFile the optional configuration file
     */
    public void setConfigFile(String configFile) {
        this.configFile = configFile;
    }

    /**
     * Returns the number of beams.
     *
     * @return the number of beams
     */
    public int getNumberOfBeams() {
        return numberOfBeams;
    }

    /**
     * Sets the number of beams.
     *
     * @param numberOfBeams the number of beams
     */
    public void setNumberOfBeams(int numberOfBeams) {
        this.numberOfBeams = numberOfBeams;
    }

    /**
     * Returns the batch size.
     *
     * @return the batch size
     */
    public int getBatchSize() {
        return batchSize;
    }

    /**
     * Sets the batch size.
     *
     * @param batchSize the batch size
     */
    public void setBatchSize(int batchSize) {
        this.batchSize = batchSize;
    }

    /**
     * Returns whether knapsack beam search is used.
     *
     * @return whether knapsack beam search is used
     */
    public boolean isUseKnapsack() {
        return useKnapsack;
    }

    /**
     * Sets whether knapsack beam search is used.
     *
     * @param useKnapsack whether knapsack beam search is used
     */
    public void setUseKnapsack(boolean useKnapsack) {
        this.useKnapsack = useKnapsack;
    }

    /**
     * Returns whether all beam search predictions are saved.
     *
     * @return whether all beam search predictions are saved
     */
    public boolean isSaveAllPredictions() {
        return saveAllPredictions;
    }

    /**
     * Sets whether all beam search predictions are saved.
     *
     * @param saveAllPredictions whether all beam search predictions are saved
     */
    public void setSaveAllPredictions(boolean saveAllPredictions) {
        this.saveAllPredictions = saveAllPredictions;
    }

    /**
     * Returns whether CPU execution is forced.
     *
     * @return whether CPU execution is forced
     */
    public boolean isForceCpu() {
        return forceCpu;
    }

    /**
     * Sets whether CPU execution is forced.
     *
     * @param forceCpu whether CPU execution is forced
     */
    public void setForceCpu(boolean forceCpu) {
        this.forceCpu = forceCpu;
    }

    /**
     * Null-safe string comparison.
     *
     * @param a the first value
     * @param b the second value
     *
     * @return true if the two values are equal
     */
    protected boolean safeEquals(String a, String b) {
        return a == null ? b == null : a.equals(b);
    }
}
