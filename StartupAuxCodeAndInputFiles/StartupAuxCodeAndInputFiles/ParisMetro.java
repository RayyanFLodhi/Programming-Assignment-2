// Paris (P2- Part B)
// Startup code given in the Fall 2025 for csi2110/csi2510
// This file only contains basic commands to read the data from the input file 'metro.txt'
// Use and modify it freely
// 
// Do not forget to add your name and student number on each program file you submit
//

import java.io.*;
import java.util.*;
import java.util.Scanner;


public class ParisMetro {

    public static void main (String[] args) {
        readMetro();
    }



    public static void readMetro()
    {

        System.out.println("\nStart Reading Metro");
        Scanner scan = new Scanner(System.in); 

        int n = scan.nextInt(); // number of vertices
        int m = scan.nextInt(); // number of edges
        System.out.println("Paris Metro Graph has "+n+" vertices and "+m+" edges.");
        

        for (int i=0; i <n; i++)
        {   int vertexNumber=Integer.valueOf(scan.next()); // vertex number (unique, example: 0016)
            String stationName=scan.nextLine(); 	    // vertex name (not unique, example: Bastille)
            //System.out.println(vertexNumber+" "+stationName);
        }

        scan.nextLine();//read the $ sign 

        for (int i=0; i <m; i++)
        {
            int v1=scan.nextInt(); int v2=scan.nextInt(); int weight=scan.nextInt(); // edge information
            //System.out.println("v1="+v1+" v2="+v2+" weight="+weight);
  
        }

       System.out.println("End Reading Metro\n");
    }


}
