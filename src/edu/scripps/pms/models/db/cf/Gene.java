package edu.scripps.pms.models.db.cf;

import java.io.Serializable;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;

/** @author Hibernate CodeGenerator */
public class Gene implements Serializable {

    /** identifier field */
    private String id;

    /** nullable persistent field */
    private String name;

    /** nullable persistent field */
    private String refseqAcc;

    /** nullable persistent field */
    private String description;

    /** full constructor */
    public Gene(String name, String refseqAcc, String description) {
        this.name = name;
        this.refseqAcc = refseqAcc;
        this.description = description;
    }

    /** default constructor */
    public Gene() {
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getName() {
        return this.name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getRefseqAcc() {
        return this.refseqAcc;
    }

    public void setRefseqAcc(String refseqAcc) {
        this.refseqAcc = refseqAcc;
    }

    public String getDescription() {
        return this.description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String toString() {
        return new ToStringBuilder(this)
            .append("id", getId())
            .toString();
    }

    public boolean equals(Object other) {
        if ( !(other instanceof Gene) ) return false;
        Gene castOther = (Gene) other;
        return new EqualsBuilder()
            .append(this.getId(), castOther.getId())
            .isEquals();
    }

    public int hashCode() {
        return new HashCodeBuilder()
            .append(getId())
            .toHashCode();
    }

}
