package qcorr;

import java.io.*;

/**
 * @author  David L. Tabb
 * @version $Id: EndsWithFilter.java,v 1.2 2007/02/27 01:31:02 taoxu Exp $
 */
/* This filter allows one to do directory listings including only
 * files with a particular filename suffix.  In this case, we want
 * to create lists of only the .out files from the subdirectories.
 */

public class EndsWithFilter implements FilenameFilter {
    private String		FilterExtension;
    public EndsWithFilter(String extension) {
	FilterExtension = extension;
    }
    public boolean accept(File dir, String name) {
	return (name.endsWith(FilterExtension));
    }
}
