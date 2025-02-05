using static System.Console;
using static System.Math;
public class vec{

public double x,y,z; //the three components of a vector

public vec(){ x=y=z=0; } //Default values are set
public vec(double x,double y,double z){
	this.x=x; this.y=y; this.z=z;
	} //Parametized constructor

public static vec operator*(vec v, double c){
	return new vec(c*v.x,c*v.y,c*v.z);
	} // omskriver gange operatoren til når en konstant ganges på en vektor

public static vec operator*(double c, vec v){
	return v*c;
	} // omkskriver gange operatoren ved at udnytte overstående definition så man både kan soge c*vec og vec*c

public static vec operator-(vec u){
	return new vec(-u.x,-u.y,-u.z);
	} // Man kan nu sætte minus foran en vektor

public static vec operator-(vec u, vec v){
	return new vec(u.x-v.x,u.y-v.y,u.z-v.z);
	} // Man trækker to vektorer fra hinanden

public static vec operator+(vec u, vec v){
	return new vec(u.x+v.x,u.y+v.y,u.z+v.z);
	} // Man kan nu lægge to vektorer til hinanden

public void print(string s=""){
	Write(s); WriteLine($" [{x}, {y}, {z}]");
	}

public double dot(vec other){
	return this.x*other.x+this.y*other.y+this.z*other.z;
}

public static double dot(vec u, vec v){
	return u.x*v.x+u.y*v.y+u.z*v.z;
}
public override string ToString(){ return $"[{x}, {y}, {z}]"; }


public static bool approx(double a,double b,double acc=1e-9,double eps=1e-9){
	if(Abs(a-b)<acc)return true;
	if(Abs(a-b)<(Abs(a)+Abs(b))*eps)return true;
	return false;
	}

public bool approx(vec other){
	if(!approx(this.x,other.x))return false;
	if(!approx(this.y,other.y))return false;
	if(!approx(this.z,other.z))return false;
	return true;
	}

public static bool approx(vec u, vec v){
	return u.approx(v);
}

}//vec
