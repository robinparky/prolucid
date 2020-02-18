
public class PhosphoResult {
    private double ascore1 = 0;
    private double ascore2 = 0;
    private double ascore3 = 0;
    private double debunkerscore = 0;
    private double confidence = 0; // dtaselect confidence score
    private double xcorr = 0;
    private double dcn = 0;
    private String sequence = "";
    private String scan = "";
    private String fragmentationMethod = "";
    private String spectrumfile = "";
    private String proteinAcc = "";
    private String proteinDescription = "";
    private String line;

    PhosphoResult(String line) {
//System.out.println("line: " + line);
        this.line = line;
        String [] arr = line.split("\t");
        if(arr.length == 13) {
            proteinAcc = arr[11];
        } else {
            proteinAcc = arr[10];
        }
    }
    public String getProteinAcc() {
        return proteinAcc;
    }
    public String getLine() {
        return line;
    }
}
