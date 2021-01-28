package com.compomics.util.experiment.identification.protein_inference.fm_index;

import com.compomics.util.experiment.io.biology.protein.Header;
import java.io.Serializable;

/**
 * Accession meta data.
 *
 * @author Dominik Kopczynski
 */
public class AccessionMetaData implements Serializable {

    private static final long serialVersionUID = 8320902054239025L;
    /**
     * The header as string.
     */
    private String headerAsString;
    /**
     * The index.
     */
    int index;
    /**
     * The index part.
     */
    int indexPart;
    /**
     * Beginning of protein sequence in the proteome.
     */
    int trueBeginning;

    /**
     * Empty default constructor.
     */
    public AccessionMetaData() {
    }
    
    /**
     * Copy constructor
     * 
     * @param accessionMetaData the item to copy
     */
    public AccessionMetaData(AccessionMetaData accessionMetaData){
        headerAsString = new String(accessionMetaData.getHeaderAsString());
        index = accessionMetaData.index;
        indexPart = accessionMetaData.indexPart;
        trueBeginning = accessionMetaData.trueBeginning;
    }

    /**
     * Constructor.
     *
     * @param header the header as parsed from the FASTA file
     */
    public AccessionMetaData(String header) {
        this.headerAsString = header;
    }

    /**
     * Constructor.
     *
     * @param header the header as parsed from the FASTA file
     * @param index the index
     * @param indexPart the index part
     * @param beginning beginning of the protein in the proteome sequence
     */
    public AccessionMetaData(String header, int index, int indexPart, int beginning) {
        this.headerAsString = header;
        this.index = index;
        this.indexPart = indexPart;
        this.trueBeginning = beginning;
    }

    /**
     * Returns the header string representation.
     *
     * @return the parsed header
     */
    public String getHeaderAsString() {
        return headerAsString;
    }
}
