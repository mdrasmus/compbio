#!/usr/bin/awk -f
#
# Author: Pouya Kheradpour 
# Mon Oct 25 15:56:09 EDT 2010
#

## INPUT FORMAT:
# Row Col Val|ColorSpace|Text|TextOpts;...
#
# for Row, Col values separated by |:
# 1)  Printed value
# 2)  Control code:
#     (empty):   standard row/column
#     Spacer:    spacer row/column. "Printed value" ignored
#     Name:      Name value centered between this and next name of same offset
#                   R/S: not rotated
#                   X: name is printed on opposite side
#                   H: put "here"
#                   B: bar is printed
#                   T: Tick mark to be printed (either "here" or center)
#     Text:      resize row/column to fit text
#                   S: name at top is rotated
#                   C/R: text is center/right on the inside
# 3)  Width (columns)/Height (rows)/Offset (Names)
# 4+) Unique identifier (used to make row/column unique)
#
# Val options (must go after number):
# -) then no box (or color)
# x) add dashed line
#
# Val can also be an (r,g,b)[(x1,x2,y1,y2)] (color in int, subbox)
#
# TextOpts
#    (r,g,b) color [default: Auto]
#     size (floating point) [default: 15 for Text; 7.5 for non-text]
#     L/C/R (left [default for text], center [default for non-text], right)
#     S/F (text rotated/not rotated)
#     NOTE: L/C/R/S only used with one Val
#
## PARAMETERS:
# C:   Colors. Space separated list of semicolon separated lists of [Val]:[Color]:(x1,x2,y1,y2)
#      Color either in Cs (below) or (r,g,b). 
#      Each value in outer list is a ColorSpace (numbered from 1)
#      (x1,x2,y1,y2) subbox interpolates for bar plot. default is 0,1,0,1
# 
# (see below for others)


END{
	OFS=FS=" "

	if (BoxHeight == "")
		BoxHeight = 20
	
	if (BoxWidth == "")
		BoxWidth = 20

	if (LFontSize == "")
		LFontSize = 15

	if (FontSize == "")
		FontSize = 15

	if (SFontSize == "")
		SFontSize = 7.5

	# replaces any missing value
	if (D == "")
		D = "x"

	if (StrokeC == "")
		StrokeC = "black"

	if (StrokeC == "-")
		StrokeC = ""

	if (SpacerHeight == "")
		SpacerHeight = 10

	if (SpacerWidth == "")
		SpacerWidth = 10

	if (LegendNum == "")
		LegendNum = 11

	if (LegendFmt == "")
		LegendFmt = "%.1f"

	# stroke around each large block
	if (StrokeSC == "-")
		StrokeSC = ""

	# auto color indicates choose black/white based on brightness
	if (STextC == "")
		STextC = "Auto"

	if (Font == "")
		Font = "Calibri"

	BoxSize[1] = BoxHeight
	BoxSize[2] = BoxWidth
	SpacerSize[1] = SpacerHeight
	SpacerSize[2] = SpacerWidth

	# specifies legend... space separated (for each color scale); _ converted to space
	LegendL = split(L, Legend, / /)

	# colors are specified val:color, where color is one of below or (r,g,b) in 0->255
	# spaces separate different color spaces
	# 4th color space from colorbrewer2.org (see Plot.py)
	if (C == "")
		C = "-1:blue;0:white;1:red -1:green;0:white;1:yellow 0:white;1:black 0:white;1:(251,128,114);2:(128,177,211);3:(179,222,105);4:(141,211,199);5:(252,205,229);6:(255,237,111);7:(217,217,217);8:(188,128,189);9:(204,235,197);10:(255,255,179);11:(253,180,98);12:(190,186,218)"

	Cs["blue"] = "(0,0,255)"
	Cs["red"] = "(255,0,0)"
	Cs["white"] = "(255,255,255)"
	Cs["yellow"] = "(255,255,0)"
	Cs["green"] = "(0,255,0)"
	Cs["black"] = "(0,0,0)"
	Cs["orange"] = "(255,127,0)"
	Cs["violet"] = Cs["purple"] = "(127,0,255)"

	if (StrokeC in Cs)
		StrokeC = "rgb" Cs[StrokeC]

	if (STextC in Cs)
		STextC = "rgb" Cs[STextC]


	################# SETUP COLOR VECTORS (INCL. SORTING THEM)
	na = split(C,A,/ /)
	for (i=1; i<=na; i++)
	{
		# insertion sort (stable)
		nC[i] = split(A[i],B,/;/)
		for (k=1; k<=nC[i]; k++)
		{
			minj = minval = 0
			for (j=k; j<=nC[i]; j++)
			{
				split(B[j], BB, /:/)

				if (minj == 0 || BB[1] < minval)
				{
					minj = j
					minval = BB[1]
				}
			}

			split(B[minj], BB, /:/)
			Cval[i, k] = BB[1]
			Ccol[i, k] = ((BB[2] in Cs) ? Cs[BB[2]] : BB[2])
			Ccol[i, k] = substr(Ccol[i, k],2,length(Ccol[i, k])-2)
			Cpos[i, k] = 3 in BB ? substr(BB[3],2,length(BB[3])-2) : "0,1,0,1"

			B[minj] = B[k]
		}
	}


	################# SETUP APPROX WIDTH OF CHARACTERS
	split("! \" # $ % & ' ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : ; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z [ \\ ] ^ _ ` a b c d e f g h i j k l m n o p q r s t u v w x y z { | }", A, / /)
	split("18 26 39 36 55 80 12 20 19 34 37 15 22 15 30 38 36 36 36 38 36 38 37 38 37 15 16 16 36 37 34 65 45 40 41 46 35 34 45 43 14 19 40 33 62 45 49 38 51 41 34 39 45 44 69 40 38 36 21 28 17 36 80 18 33 38 31 36 36 25 36 37 13 14 35 13 58 37 39 38 36 27 28 25 36 35 55 33 35 29 23 21 24", B, / /)
	for (i in A)
		CharWidth[A[i]] = B[i]
	CharWidth[" "] = 12


	# ZP[1,idx] --> Y coordinate of top most position of row idx
	# ZP[2,idx] --> X coordinate of left most position of column idx
	for (d=1; d<=2; d++)
	{
		################## GET OFFSETS FOR ROWS/COLUMNS
		# idx==0 is offset before row/column headers
		for (z=1; z<=NumZ[3-d]; z++) 
			if (ZEx[3-d,z] ~ /^Name/ && !H("X", ZEx[3-d,z], 4))
				ZP[d,0] = max(ZP[d,0], 10 + (H("RS", ZEx[3-d,z],4) ? strwidth(ZN[3-d,z]) : 20))

		# idx==1 is offset before table itself
		for (z=1; z<=NumZ[3-d]; z++) 
			if ((ZEx[3-d,z] == "" || (ZEx[3-d,z] ~ /^Text/ && !H("S", ZEx[3-d,z],4))) && ZN[3-d,z] != "")
				ZP[d,1] = max(ZP[d,1], strwidth(ZN[3-d,z]) + 20)

		ZP[d,1] += ZP[d,0]
		ZP[d,1] = max(ZP[d,1], 20)

		################### COMPUTE POSITIONS OF ROWS/COLUMNS
		for (z=1; z<=NumZ[d]; z++)
		{
			ZP[d,z+1] = ZP[d,z]

			if (ZS[d,z] != "" && ZEx[d,z]  !~ /^Name/)
				ZP[d,z+1] += ZS[d,z]
			else if (ZEx[d,z] == "")
				ZP[d,z+1] += BoxSize[d]
			else if (ZEx[d,z] == "Spacer")
				ZP[d,z+1] += SpacerSize[d]
			else if (ZEx[d,z] ~ /^Text/)
			{
				w = 0
				for (zx=1; zx<=NumZ[3-d]; zx++) 
					if (ZEx[3-d,zx] == "" || ZEx[3-d,zx] ~ /^Text/)
					{
						split(d==1 ? Val[z,zx] : Val[zx,z], A, /[|]/)
						w = max(w, strwidth(A[3]) + 10)
					}
				ZP[d,z+1] += w
			}
		}

		# idx==NumZ[d]+1 is end of table; idx==NumZ[d]+2 is end of table + NameXs
		for (z=1; z<=NumZ[3-d]; z++) 
			if (ZEx[3-d,z] ~ /^Name/ && H("X", ZEx[3-d,z], 4))
				ZP[d,NumZ[d]+2] = max(ZP[d,NumZ[d]+2], 10 + (H("RS", ZEx[3-d,z],4) ? strwidth(ZN[3-d,z]) : 10))
	
		ZP[d,NumZ[d]+2] += ZP[d,NumZ[d]+1]
	}

	LegendY[0] = ZP[1,NumZ[1]+2]
	for (i=1; i<=LegendL; i++)
	{
		LegendY[i] = LegendY[i-1]
		if (Legend[i] != "")
			LegendY[i] += 3 * 20
	}

	################## PRINT SVG HEADER (includes size)
	print "<?xml version=\"1.0\" standalone=\"no\"?><!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
	print "<svg width=" Q(20 + max(ZP[2,NumZ[2]+2], ZP[2,1] + ((L == "") ? 0 : (LegendNum * 20)))) " height=" Q(20 + LegendY[LegendL]) " version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">"

	################# PRINT NAMES
	for (d=1; d<=2; d++)
		for (z=1; z<=NumZ[d]; z++)
			if (ZEx[d,z] ~ /^Name/)
			{
				isX = H("X", ZEx[d,z], 4)
				isR = H("RS", ZEx[d,z], 4)
				isH = H("H", ZEx[d,z], 4)

				for (zx=z+1; zx<=NumZ[d]; zx++)
					if (ZEx[d,zx] ~ /^Name/ && isX == H("X", ZEx[d,zx], 4) && ZS[d,zx] == ZS[d,z])
						break

				p = isH ? (ZP[d,z]+4) : (ZP[d,z] + ZP[d,zx])/2
				m = isX ? -1 : 1
				px = ZS[d,z] + ZP[3-d,isX ? (NumZ[3-d]+1) : 0] - 17 * m
				o = isH ? 1 : (isR ? (isX ? "" : 1) : 2)

				print text(d, ZN[d,z], px, p, isR ? "" : 1, o, LFontSize)

				# print bar
				if (H("B", ZEx[d,z], 4))
					print poly(d, (px + LFontSize*m*2/3) " " (px + LFontSize*m*2/3), (ZP[d,z] + 3) " " (ZP[d,zx] - 3), "", "black")

				# print tick mark
				if (H("T", ZEx[d,z], 4))
				{
					if (isH)
						print poly(d, (px - 10) " " (px + 10), ZP[d,z] " " ZP[d,z], "", "black")
					else
						print poly(d, (px + 6*m) " " (px + 14*m), (p " " p), "", "black")
				}
			}

	################## PRINT HEADERS
	for (d=1; d<=2; d++)
		for (z=1; z<=NumZ[d]; z++)
			if (ZEx[d,z] == "" || (ZEx[d,z] ~ /^Text/ && !H("S", ZEx[d,z], 4)))
				print text(d, ZN[d,z], ZP[3-d,1] - 4, (ZP[d,z] + ZP[d,z+1])/2, "", 1)
			else if (ZEx[d,z] ~ /^Text/)
				print text(d, ZN[d,z], ZP[3-d,1] - 11, ZP[d,z] + 3, 1, 1)


	################# PRINT MATRIX
	for (r=1; r<=NumZ[1]; r++)
		for (c=1; c<=NumZ[2]; c++)
			if ((ZEx[1,r] == "" || ZEx[1,r] ~ /^Text/) && (ZEx[2,c] == "" || ZEx[2,c]  ~ /^Text/))
				print rect((Val[r,c] == "") ? D : Val[r,c], ZP[2,c], ZP[1,r], ZP[2,c+1], ZP[1,r+1], 
					(ZEx[1,r] ~ /^Text/ || ZEx[2,c] ~ /^Text/) ? FontSize : SFontSize,
					((ZEx[1,r] ~ /^Text/ || ZEx[2,c] ~ /^Text/) ? "L" : "C") ((ZEx[1,r]  ~ /^Text/) ? "S" : "F"))


	################# PRODUCE STROKES AROUND LARGE BLOCKS
	if (StrokeSC != "")
	{
		for (r1=1; r1<NumZ[1]; r1++)
		{
			# find a spacer
			for (r2=r1+1; r2<=NumZ[1] && (ZEx[1,r2]  != "Spacer"); r2++);

			for (c1=1; c1<=NumZ[2]; c1++)
			{
				for (c2=c1+1; c2<=NumZ[2] && (ZEx[2,c2]  != "Spacer"); c2++);
				print poly(1, ZP[2,c1]" "ZP[2,c2]" "ZP[2,c2]" "ZP[2,c1], ZP[1,r1]" "ZP[1,r1]" "ZP[1,r2]" "ZP[1,r2], "", StrokeSC)
				c1 = c2
			}
			r1 = r2
		}
	}


	################# PRINT LEGEND
	# always have stroke color black
	StrokeC = "black"
	for (i=1; i<=LegendL; i++)
		if (Legend[i] != "")
		{
			print text(1, gensub(/_/, " ", "g", Legend[i]), ZP[2,1] + 20 * LegendNum/2, LegendY[i]-3*20/2, "", 2)
			for (j=1; j<=LegendNum; j++)
				print rect((Cval[i, 1] + (j-1)/(LegendNum-1) * (Cval[i, nC[i]] - Cval[i, 1])) "|" i "|" sprintf(LegendFmt, (Cval[i, 1] + (j-1)/(LegendNum-1) * (Cval[i, nC[i]] - Cval[i, 1]))), ZP[2,1] + 20*(j-1), LegendY[i]-20, ZP[2,1] + 20*j, LegendY[i])
		}


	################# END OF SVG
	print "</svg>"
}

function color(val, sc,   i, A, B, w) 
{
	if (sc == "")
		sc = 1

	# if in (r,g,b) style, just return (stripping any extra characters)
	if (val ~ /^\([0-9]*,[0-9]*,[0-9]*\)/)
		return "rgb" gensub(/\).*/, ")", 1, val)

	if (val ~ /^[-x]*$/ || !(sc in nC))
		return "" 

	for (i=1; i<=nC[sc] && (val > Cval[sc, i]); i++);
	
	split(Ccol[sc, i-1], A, /,/)
	split(Ccol[sc, i], B, /,/)

	w = (i == 1) ? 1 : ((i > nC[sc]) ? 0 : ((val - Cval[sc, i-1]) / (Cval[sc, i] - Cval[sc, i-1])))

	return "rgb(" int(A[1] * (1-w) + B[1] * w) "," int(A[2] * (1-w) + B[2] * w) "," int(A[3] * (1-w) + B[3] * w) ")"
}

# very similar to color, but interpolates the positions
function subbox(val, sc,   i, A, B, w) 
{
	if (sc == "")
		sc = 1

	# if in (r,g,b)(x1,x2,x3,x4) style, just return bounding box
	if (val ~ /^\([0-9]*,[0-9]*,[0-9]*\)\([0-9.]*,[0-9.]*,[0-9.]*,[0-9.]*\)[-x]*$/)
		return gensub(/\).*/, "", 1, gensub(/^\([0-9]*,[0-9]*,[0-9]*\)\(/, "", 1, val))

	# if blank, or rgb style without subbox or invalid sc, then return entire box
	if (val ~ /^[-x]*$/ || val ~ /^\([0-9]*,[0-9]*,[0-9]*\)[-x]*$/ || !(sc in nC))
		return "0,1,0,1"

	for (i=1; i<=nC[sc] && (val > Cval[sc, i]); i++);
	w = (i == 1) ? 1 : ((i > nC[sc]) ? 0 : ((val - Cval[sc, i-1]) / (Cval[sc, i] - Cval[sc, i-1])))

	split(Cpos[sc, i-1], A, /,/)
	split(Cpos[sc, i], B, /,/)

	return (A[1] * (1-w) + B[1] * w) "," (A[2] * (1-w) + B[2] * w) "," (A[3] * (1-w) + B[3] * w) "," (A[4] * (1-w) + B[4] * w)
}

# x2, y2 are optional
function rect(val, x1, y1,   x2, y2, defTS, defTO,    A, n, xm, ym, B, C, P, i, T, L, NoStr, Cross, TC, TS, nn, TO, tx, ty, to, ta, x1p, x2p, y1p, y2p) 
{
	n = max(split(val,A,/[;]/),1)

	if (x2 == "") x2 = x1 + BoxSize[2]
	if (y2 == "") y2 = y1 + BoxSize[1]
	
	if (defTS == "") defTS = SFontSize
	if (defTO == "") defTO = "C"

	xm = (x1 + x2)/2
	ym = (y1 + y2)/2

	for (i=1; i<=n; i++) 
	{
		nn = split(A[i],B,/[|]/)

		C[i] = color(B[1], B[2])
		P[i] = subbox(B[1], B[2])
		TC[i] = tcolor(C[i])
		TS[i] = defTS
		TO[i] = defTO

		if (nn>=3) 
		{
			T[i] = B[3]

			if (B[4] ~ /^\([0-9]*,[0-9]*,[0-9]*\)/)
				TC[i] = "rgb" gensub(/\).*/, ")", 1, B[4])

			if (gensub(/^\([0-9]*,[0-9]*,[0-9]*\)/, "", 1, B[4]) ~ /^[0-9.]/)
				TS[i] = gensub(/^\([0-9]*,[0-9]*,[0-9]*\)([0-9.]*)[A-Za-z]*/, "\\1", 1, B[4]) + 0

			TO[i] = TO[i] gensub(/^\([0-9]*,[0-9]*,[0-9]*\)[0-9.]*/, "", 1, B[4])
		}

		if (B[1] ~ /-[x]*$/)
			NoStr = 1

		if (B[1] ~ /x[-]*$/)
			Cross = 1
	}

	if (n==1 && C[1] != "")
	{
		# only do interpolating for this one (for now)
		split(P[1], A, /,/)
		
		x1p = (x2 - x1) * A[1] + x1
		x2p = (x2 - x1) * A[2] + x1
		y1p = (y2 - y1) * A[3] + y1
		y2p = (y2 - y1) * A[4] + y1

		L = L poly(1, x1p" "x2p" "x2p" "x1p, y1p" "y1p" "y2p" "y2p, C[1])
	}
	else if (n>1)
	{
		if (n==2)
			L = L poly(1, x1" "x2" "x1, y1" "y1" "y2, C[1]) poly(1, x2" "x1" "x2, y1" "y2" "y2, C[2]) 
		else if (n==3)
			L = L poly(1, x1" "xm" "x2, y1" "y2" "y1, C[1]) poly(1, x1" "x1" "xm, y1" "y2" "y2, C[2]) poly(1, x2" "x2" "xm, y1" "y2" "y2, C[3])
		else
			L = L poly(1, x1" "x2" "xm, y1" "y1" "ym, C[1]) poly(1, x1" "x1" "xm, y1" "y2" "ym, C[2]) poly(1, x2" "x2" "xm, y1" "y2" "ym, C[3]) poly(1, x1" "x2" "xm, y2" "y2" "ym, C[4])
	}

	if (!NoStr)
		L = L poly(1, x1" "x2" "x2" "x1, y1" "y1" "y2" "y2, "", StrokeC)

	if (Cross)
		L = L poly(1, x1 " " x2, y2 " " y1, "", TC[1])

	# produce text
	if ((1 in T) && (2 in T) && (3 in T) && (4 in T))
	{
		L = L text(1, T[1], xm, y1 * 0.75 + y2 * 0.25, "", 2, TS[1], TC[1])
		L = L text(1, T[2], x1 * 0.75 + x2 * 0.25, ym, "", 2, TS[2], TC[2])
		L = L text(1, T[3], x1 * 0.25 + x2 * 0.75, ym, "", 2, TS[3], TC[3])
		L = L text(1, T[4], xm, y1 * 0.25 + y2 * 0.75, "", 2, TS[4], TC[4])
	}
	else if ((1 in T) && (2 in T) && (3 in T))
	{
		L = L text(1, T[1], xm, y1 * 0.65 + y2 * 0.35, "", 2, TS[1], TC[1])
		L = L text(1, T[2], x1 * 0.8 + x2 * 0.2, y1 * 0.25 + y2 * 0.75, "", 2, TS[2], TC[2])
		L = L text(1, T[3], x1 * 0.2 + x2 * 0.8, y1 * 0.25 + y2 * 0.75, "", 2, TS[3], TC[3])
	}
	else if ((1 in T) && (2 in T))
	{
		L = L text(1, T[1], x1 * 0.75 + x2 * 0.25, y1 * 0.75 + y2 * 0.25, "", 2, TS[1], TC[1])
		L = L text(1, T[2], x1 * 0.25 + x2 * 0.75, y1 * 0.25 + y2 * 0.75, "", 2, TS[2], TC[2])
	}
	else
	{
		# text orientation (1 -- vertical, "" -- horizontal)
		to = TO[1] !~ /S[^F]*$/ ? "" : 1

		# text alignment ("" -- left, 1 -- right; 2 -- center)
		ta = TO[1] ~ /C[^LR]*$/ ? 2 : (TO[1] ~ /R[^LC]*$/ ? 1 : "")

		# text x and y
		tx = (to || ta==2) ?  xm : (ta ? (x2-3) : (x1+3))
		ty = (!to || ta==2) ? ym : (ta ? (y1+3) : (y2-3))

		L = L text(1, T[1], tx, ty, to, ta, TS[1], TC[1])
	}
	return L
}

function tcolor (bgc,  A)
{
	if (STextC != "Auto")
		return STextC
	
	if (substr(bgc,1,3) != "rgb")
		return "rgb" Cs["black"]

	split(substr(bgc,5,length(bgc)-5), A, /,/)

	# brightness formula found on web
	return "rgb" Cs[((A[1] * 299 + A[2] * 587 + A[3] * 114) < (127500)) ? "white" : "black"]
}

function strwidth(x,  i, l)
{
	for (i=1; i<=length(x); i++)
		l += CharWidth[substr(x,i,1)]
	return l * FontSize/50
}

########################## SIMPLE HELPER FUNCTIONS
function H(x, str, skip) { return substr(str, skip+1) ~ ("[" x "]") }
function Q(x) { return "\"" x "\"" }
function max(a,b) { return a > b ? a : b }
function X(x,y,def) { return x != "" ? y : def }

########################## FUNCTIONS THAT ACTUALLY WRITE SVG
# align "" --> left
# align 1  --> right
# align 2  --> center
function text (dim, t, x, y, rot, align, size, color,  temp) 
{
	if (dim == 2)
	{
		temp = x
		x = y
		y = temp
		rot = rot ? "" : 1
		align = align == 2 ? align : (align == "" ? 1 : "")
	}

	if (size == "") size = FontSize
	if (rot)
		x += (size+1)/3
	else
		y += (size+1)/3

	gsub(/</, "\\&lt;", t)
	gsub(/>/, "\\&gt;", t)

	return X(t,"<text x=" Q(x) " y=" Q(y) " font-family=" Q(Font) " font-size=" Q(size) X(color, " fill=" Q(color)) X(rot, " transform=" Q("rotate(270 " x " " y ")")) X(align, " text-anchor=" Q(align==1 ? "end" : (align==2 ? "middle" : "start"))) ">" t "</text>")
}


function poly (dim,  xs, ys, cl, str,  X, Y, n, i, L)
{
	n = split(dim==2 ? ys : xs, X, / /)
	    split(dim==2 ? xs : ys, Y, / /)

	for (i=1; i<=n; i++)
		L = L " " X[i] "," Y[i]
	
	return X(cl str, "<polygon points=" Q(substr(L,2)) " fill=" Q(X(cl,cl,"none")) X(str, " stroke=" Q(str) " stroke-width=" Q(1)) "/>")
}


########################## READ IN STDIN
{
	for (d=1; d<=2; d++)
	{
		if ($d == "")
			$d = "|-"

		if (!((d,$d) in Zinv))
		{
			split($d,A,/[|]/)
			ZN[d,++NumZ[d]] = A[1]
			ZEx[d,NumZ[d]] = A[2]
			ZS[d,NumZ[d]] = A[3]
			Zinv[d,$d] = NumZ[d]
		}
	}

	Val[Zinv[1,$1],Zinv[2,$2]] = $3
}

