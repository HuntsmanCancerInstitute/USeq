package trans.cel;
import java.io.*;

/**Class holding info relating to one quantile value. Stripped for size.*/
public class MultiSetQuantile implements Serializable{
	
	//fields
	float value = 0;
	byte chromosomeNumber = 0;
	int position = 0;
	
	//constructor
	public MultiSetQuantile(byte chromosomeNumber, int position, float value){
		this.chromosomeNumber = chromosomeNumber;
		this.position = position;
		this.value = value;
	}
	
	public String toString(){
		return chromosomeNumber+"\t"+position+"\t"+value;
	}
}
