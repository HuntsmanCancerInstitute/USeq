package trans.cel;

/**Class holding info relating to one quantile value.*/
public class Quantile {
	
	//fields
	public float value = 0;
	public int position = 0;
	
	//constructor
	public Quantile(int position, float value){
		this.position = position;
		this.value = value;
	}
}
