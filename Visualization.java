// Software: Visualization for Genetic Algorithm to solve Travelling Saleman Problem (TSP)
// Author: Hy Truong Son
// Major: BSc. Computer Science
// Class: 2013 - 2016
// Institution: Eotvos Lorand University
// Email: sonpascal93@gmail.com
// Website: http://people.inf.elte.hu/hytruongson/
// Final update: October 4th, 2015
// Copyright 2015 (c) Hy Truong Son. All rights reserved. Only use for academic purposes.

import java.lang.*;
import java.io.*;
import java.awt.*;
import javax.swing.*;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Scanner;

import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;

public class Visualization {
	
	static int Max_nPoints = 10000;
	static int Circle_Radius = 6;
	static int scaled_ratio = 8;
	static int margin = 50;
	
	static String ImageName = "Solution.png";
	
	static int width, height;
	static String DataName, SolutionName;
	static int nPoints;
	static class aPoint {
		public int x, y;
	}
	static aPoint Point[] = new aPoint [Max_nPoints];
	static int Path[];
	
	public static String GetType(String FileName){
		int i, j;
		String res;
		j = 0;
		for (i = 0; i < FileName.length(); i++)
			if (FileName.charAt(i) == '.'){
				j = i;
				break;
			}
		res = "";
		for (i = j + 1; i < FileName.length(); i++) res += FileName.charAt(i);
		return res;
	}
	
	public static void ReadData(String FileName) throws IOException {
		BufferedReader file = new BufferedReader(new FileReader(FileName));
		
		String s, temp;
		nPoints = 0;
		
		while (true){
			s = file.readLine();
			if (s.charAt(0) == 'E') break;
			nPoints++;
			
			int t = 0;
			int i = 0;
			int count = 0;
			
			while (i < s.length()){
				if (s.charAt(i) == ' '){
					i++;
					continue;
				}
				
				temp = "";
				for (int j = i; j < s.length(); j++)
					if (s.charAt(j) != ' '){
						i = j;
						temp += s.charAt(j);
					}else break;
					
				i++;
				count++;
				
				if (count == 1){
					t = Integer.parseInt(temp);
					Point[t] = new aPoint();
				}else
					if (count == 2) Point[t].x = Integer.parseInt(temp); else
						Point[t].y = Integer.parseInt(temp);
			}
		}
		
		file.close();
	}
	
	public static void ReadSolution(String FileName) throws IOException {
		BufferedReader file = new BufferedReader(new FileReader(FileName));
		
		Path = new int [nPoints];
		file.readLine();
		for (int i = 0; i < nPoints; i++)
			Path[i] = Integer.parseInt(file.readLine());
		
		file.close();
	}
	
	public static void Draw_Circle(BufferedImage image, int x, int y, int r){
		Graphics2D g2d = image.createGraphics();
		g2d.setColor(Color.GREEN);
		y = height - y - 1;
		x -= r / 2;
		y -= r / 2;
		Ellipse2D.Double circle = new Ellipse2D.Double(x, y, r, r);
		g2d.fill(circle);
	}
	
	public static void Draw_Line(BufferedImage image, int x1, int y1, int x2, int y2){
		Graphics2D g2d = image.createGraphics();
		g2d.setColor(Color.RED);
		y1 = height - y1 - 1;
		y2 = height - y2 - 1;
		g2d.drawLine(x1, y1, x2, y2);
	}
	
	public static int RGB(int red,int green,int blue){
		return (0xff000000) | (red << 16) | (green << 8) | blue;
	}
	
	public static void Visual(String FileName) throws IOException {
		width = 0;
		height = 0;
		for (int i = 1; i <= nPoints; i++){
			width = Math.max(width, scaled_ratio * Point[i].x + 2 * margin);
			height = Math.max(height, scaled_ratio * Point[i].y + 2 * margin);
		}
	
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
				image.setRGB(i, j, RGB(255, 255, 255));
				
		for (int i = 0; i < nPoints; i++){
			int x1 = margin + scaled_ratio * Point[Path[i]].x;
			int y1 = margin + scaled_ratio * Point[Path[i]].y;
			int x2 = margin + scaled_ratio * Point[Path[(i + 1) % nPoints]].x;
			int y2 = margin + scaled_ratio * Point[Path[(i + 1) % nPoints]].y;
			
			Draw_Line(image, x1, y1, x2, y2);
		}
		
		for (int i = 1; i <= nPoints; i++)
			Draw_Circle(image, margin + scaled_ratio * Point[i].x, margin + scaled_ratio * Point[i].y, Circle_Radius);
		
		File file = new File(FileName);
		ImageIO.write(image, GetType(FileName), file);
	}
	
	public static void main(String args[]) throws IOException{
		BufferedReader Buffer = new BufferedReader(new InputStreamReader(System.in));
		System.out.print("Data file contains points' coordinates: ");
		DataName = Buffer.readLine();
		System.out.print("Solution file: ");
		SolutionName = Buffer.readLine();
		
		ReadData(DataName);
		ReadSolution(SolutionName);
		Visual(ImageName);
	}	
	
}
