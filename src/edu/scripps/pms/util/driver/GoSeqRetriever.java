/*
 * @(#)GoSeqRetriever.java
 *
 * Copyright Notice:
 *
 * Copyright 2004 UBC Bioinformatics Centre
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 *        Free Software Foundation, Inc.,
 *        59 Temple Place, Suite 330, Boston, MA
 *        02111-1307  USA
 * or visit http://www.gnu.org/copyleft/gpl.html
 *
 */

/**
 * @file GoSeqRetriever.java
 * This is the source file for ubic.atlas.locuslink.db.GoSeqRetriever
 *
 * @author Tao Xu
 * @date $Date
 */


import edu.scripps.pms.go.GoGet;






/**
 * @class GoSeqRetriever
 * This class is a subclass of MysqlDb. It provides an interface to
 * access GO database  
 *
 */
public class GoSeqRetriever {
    public static final String USAGE = "!!! USAGE: java -Xmx1000M dbname outputFileName !!!";
    public static void main(String [] args) throws Exception {
        //String pid = args[0];
        try { 
            GoGet gg = new GoGet(args[0]);
            gg.outputAllProteinSeqInFasta(args[1]);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
}





