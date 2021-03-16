package com.compomics.util.experiment.io.biology.protein;

/**
 * Interface for a class able to provide metadata on proteins.
 *
 * @author Marc Vaudel
 */
public interface ProteinDetailsProvider {

    /**
     * Returns the description of the protein with the given accession.
     *
     * @param accession the accession of the protein
     *
     * @return the description of the protein with the given accession
     */
    public String getDescription(String accession);

    /**
     * Returns the simple description of the protein with the given accession.
     *
     * @param accession the accession of the protein
     *
     * @return the description of the protein with the given accession
     */
    public String getSimpleDescription(String accession);

    /**
     * Returns the the protein database for the given protein.
     *
     * @param accession the accession of the protein
     *
     * @return the name of the protein database
     */
    public ProteinDatabase getProteinDatabase(String accession);

    /**
     * Returns the gene name for the given protein.
     *
     * @param accession the accession of the protein
     *
     * @return the gene name for the given protein
     */
    public String getGeneName(String accession);

    /**
     * Returns the taxonomy for the given protein.
     *
     * @param accession the accession of the protein
     *
     * @return the taxonomy for the given protein
     */
    public String getTaxonomy(String accession);
    
    /**
     * Returns the organism identifier for the given protein.
     *
     * @param accession the accession of the protein
     *
     * @return the organism name for the given protein
     */
    public String getOrganismIdentifier(String accession);

    /**
     * Returns an integer representing the protein evidence level as indexed by
     * UniProt.
     *
     * @param accession the protein accession
     *
     * @return an integer representing the protein evidence level as indexed by
     * UniProt
     */
    public Integer getProteinEvidence(String accession);

}
