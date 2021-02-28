/*
 * 天文計算関係スクリプト version 0.18j at 2021/02/27
 * Copyright (c) 1999-2001, 2004, 2005, 2017, 2021 Yoshihiro Sakai & Sakai Institute of Astrology
 * This software is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 * 2017/06/05[0.16j] 軌道要素６パラ版に対応
 * 2017/09/15[0.17j] ΔＴの計算式を見直すついでに出典を書く
 * 2021/02/27[0.18j] 二体問題まわりロジック見直し
 */

// グレゴリオ暦専用！
function cnvCalendar( JD ){
	JD += 0.5;
	var Z = Math.floor( JD );
	var F = JD - Z;

	var A = 0;
	if(Z >= 2299161){
		var alpha = Math.floor((Z - 1867216.25) / 36524.25);
		A = Z + 1 + alpha - Math.floor(alpha / 4);
	} else {
		A = Z;
	}

	var B = A + 1524;
	var C = Math.floor((B - 122.1) / 365.25);
	var D = Math.floor(365.25 * C);
	var E = Math.floor((B - D) / 30.6001);

	var da = B - D - Math.floor(30.6001 * E);
	var mo = (E < 13.5) ? (E - 1) : (E - 13);
	var ye = (mo > 2.5) ? (C - 4716) : (C - 4715);

	var ti = F * 24.0;
	var ho = Math.floor(ti);
	var mi = (ti - ho) * 60.0;

	var res = new Array(ye, mo, da, ho, mi);
	return res;
}

// その日のユリウス日を計算する
function calJD(ye, mo, da, ho, mi){ // 実数体上
	var y0 = (mo > 2) ? ye : (ye -  1);
	var m0 = (mo > 2) ? mo : (mo + 12);
	var JD = Math.floor(365.25 * y0) + Math.floor(y0 / 400) - Math.floor(y0 / 100);
	JD	+= Math.floor(30.59 * (m0 - 2)) + da;
	JD	+= ((ho - 9) * 60.0 + mi) / 1440.0 + 1721088.5;

	return JD;
}

function cnvJDr( JD ){ // 実数体→整数環
	var date = cnvCalendar(JD);
	var ye = date[ 0 ];
	var mo = date[ 1 ];
	var da = date[ 2 ];
	var JDz = calJDz(ye, mo, da);

	return JDz;
}

function calJDz(year, month, day){ // 整数環上
	var yt = year;
	var mt = month;
	var dt = day;

	if(month < 3){
		yt--;
		mt += 12;
	}

	var JD = Math.floor(365.25 * yt) + Math.floor(30.6001 * (mt + 1));
	JD += dt + 1720995;
	JD += 2 - Math.floor(yt / 100) + Math.floor(yt / 400);

	return JD;
}

// d, Tを配列で返す。
function calTimeCoefficient( JD ){
	var d = JD - 2451545.0;
	var T = d / 36525.0;
	var coef = new Array(d, T);

	return coef;
}

// 軌道要素６パラ版→通常の軌道要素
// a : semi-major axis (au).
// l : mean longitude (degree).
// h : e * sin(pi).
// k : e * cos(pi).
// p : g * sin(om).
// q : g * cos(om).
// e : eccentricity.
// g : sine of the half inclination.
// pi: longitude of the perihelion.
// om: longitude of the ascending node.
function convertOrbitalElement( a, l, h, k, p, q ) {
	var L = l;
	var e = Math.sqrt( h * h + k * k );
	var g = Math.sqrt( p * p + q * q );
	var i = Math.asin( g ) * 2.0 / deg2rad;
	var opi = Math.atan2( h, k ) / deg2rad;
	var omg = Math.atan2( p, q ) / deg2rad;

	var result = [L, opi, omg, i, e, a];
	return result;
}

// まとめて軌道計算。返値は（黄経、黄緯、動径）。
function orbitWork(L, opi, omg, i, e, a){
	var M = mod360(L - opi);
	var E = mod360(solveKepler(M, e));

	var rE  = E * deg2rad;
	var thv = Math.sqrt( (1 + e) / (1 - e) ) * Math.tan(rE / 2.0);
	var v   = mod360(Math.atan2(thv, 1.0) * 2.0 / deg2rad);
	var r   = a * (1.0 - e * Math.cos(rE));
	var u   = L + v - M - omg;
	
	var ri = i * deg2rad;
	var ru = u * deg2rad;
	var l  = mod360(omg + Math.atan2(Math.cos(ri) * Math.sin(ru), Math.cos(ru)) / deg2rad);

	var rb = Math.asin(Math.sin(ru) * Math.sin(ri));
	var b  = rb / deg2rad;

	var res = new Array(l, b, r);
	return res;
}

// Kepler方程式(M = E - e sinE)を解く。
function solveKepler(M, e){
	var Mr = M * deg2rad;
	var Er = Mr;

	for( var i = 0; i < 20; i++ ){
		Er = Mr + e * Math.sin( Er );
	}

	return Er / deg2rad;
}

// 日心位置から地心位置へコンバートし、地心黄経を返す。
function convertGeocentric( earthCoor, planetCoor ){
	var lp, bp, rp, ls, bs, rs;

	ls = earthCoor[ 0 ];
	bs = earthCoor[ 1 ];
	rs = earthCoor[ 2 ];

	lp = planetCoor[ 0 ];
	bp = planetCoor[ 1 ];
	rp = planetCoor[ 2 ];

	var xs = rs * Math.cos(ls * deg2rad) * Math.cos(bs * deg2rad);
	var ys = rs * Math.sin(ls * deg2rad) * Math.cos(bs * deg2rad);
	var zs = rs *                          Math.sin(bs * deg2rad);
	var xp = rp * Math.cos(lp * deg2rad) * Math.cos(bp * deg2rad);
	var yp = rp * Math.sin(lp * deg2rad) * Math.cos(bp * deg2rad);
	var zp = rp *                          Math.sin(bp * deg2rad);

	var xg = xp + xs;
	var yg = yp + ys;
	var zg = zp + zs;

	var rg = Math.sqrt(xg * xg + yg * yg + zg * zg);
	var lg = Math.atan2(yg, xg) / deg2rad;
	if(lg < 0.0) lg += 360.0;
	var bg = asin4deg(zg / rg);

	var res = new Array(lg, bg, rg);
	return res;
}

// 黄道座標系から赤道座標系へ変換する。
function convertEquatorial(lon, lat, obl){
	var xs = cos4deg(lon) * cos4deg(lat);
	var ys = sin4deg(lon) * cos4deg(lat);
	var zs =                sin4deg(lat);

	var xd = xs;
	var yd = ys * cos4deg(obl) - zs * sin4deg(obl);
	var zd = ys * sin4deg(obl) + zs * cos4deg(obl);

	var RA = atan2(yd, xd) / deg2rad;
	if(RA < 0.0) RA += 360.0;
	var Dec = asin4deg(zd);

	var res = new Array(RA, Dec);
	return res;
}

// 歳差補正
function coordinateConvertFromJ2000( arg ){
	var zeta, zz, theta;
	var x, y, z, xd, yd, zd, xs, ys, zs;

	var xs  = arg[ 0 ];
	var ys  = arg[ 1 ];
	var zs  = arg[ 2 ];
	var tjd = arg[ 3 ];

	var T = (tjd - 2451545.0) / 36525.0;

	var zeta  = ((( 0.017998 * T + 0.30188) * T + 2306.2181) * T) / 3600.0;
	var zz    = ((( 0.018203 * T + 1.09468) * T + 2306.2181) * T) / 3600.0;
	var theta = (((-0.041833 * T - 0.42665) * T + 2004.3109) * T) / 3600.0;

// Step 1
	x =  sin4deg(zeta) * xs + cos4deg(zeta) * ys;
	y = -cos4deg(zeta) * xs + sin4deg(zeta) * ys;
	z = zs;

// Step 2;
// 	x = x;
	y = 0 * x + cos4deg(theta) * y + sin4deg(theta) * z;
	z = 0 * x - sin4deg(theta) * y + cos4deg(theta) * z;

// Step 3
	xd = -sin4deg(zz) * x - cos4deg(zz) * y;
	yd =  cos4deg(zz) * x - sin4deg(zz) * y;
	zd = z;

	var res = new Array(xd, yd, zd);
	return res;
}

// 地方恒星時計算
function calLST(JD, ho, mi, lo){
	var JD0 = Math.floor(JD - 0.5) + 0.5;
	var T = (JD0 - 2451545.0) / 36525.0;
	var UT = (JD - JD0) * 360.0 * 1.002737909350795;
	if( UT < 0 ) UT += 360.0;

	//グリニッジ恒星時計算
	var GST  = 0.279057273 + 100.0021390378 * T + 1.077591667e-06 * T * T;
	    GST  = GST - Math.floor(GST);
	    GST *= 360.0;

	// 地方恒星時計算＋章動補正
	var LST = mod360(GST + UT + lo);
	var dpsi = calNutation(JD);
	var eps  = calOblique(JD);
	LST += dpsi * cos4deg(eps) / 3600.0;
	if(LST < 0.0) LST += 360.0;

	return LST;
}

// 黄道傾斜角を計算する関数
function calOblique(JD ){
	var T = (JD - 2451545.0) / 36525.0;
	var Omg = mod360(125.00452 - T *   1934.136261);
	var Ls  = mod360(280.4665  + T *  36000.7698);
	var Lm  = mod360(218.3165  + T * 481267.8813);

	var e = 84381.448 + T * (-46.8150 + T * (-0.00059 + T * 0.001813));
	var deps  =  9.20 * cos4deg(1.0 * Omg);
	    deps +=  0.57 * cos4deg(2.0 * Ls);
	    deps +=  0.10 * cos4deg(2.0 * Lm);
	    deps += -0.09 * cos4deg(2.0 * Omg);

	return (e + deps) / 3600.0;
}

// 章動を計算する関数（簡略版）
function calNutation( JD ){
	var T = (JD - 2451545.0) / 36525.0;

	var Omg = mod360(125.00452 - T *   1934.136261);
	var Ls  = mod360(280.4665  + T *  36000.7698);
	var Lm  = mod360(218.3165  + T * 481267.8813);

	var dpsi  = -17.20 * sin4deg(1.0 * Omg);
	    dpsi +=  -1.32 * sin4deg(2.0 * Ls);
	    dpsi +=  -0.23 * sin4deg(2.0 * Lm);
	    dpsi +=   0.21 * sin4deg(2.0 * Omg);

	return dpsi;
}

// 均時差を計算する関数
function calEqT( JD ){
	var  T = ( JD - 2451545.0 ) / 36525.0;

	var L0  = mod360( 36000.76983 * T );
		L0  = mod360( 280.46646 + L0 + 0.0003032 * T * T );
		L0 *= deg2rad;

	var  M  = mod360( 35999.05029 * T );
		 M  = mod360( 357.52911 + M  - 0.0001537 * T * T );
		 M *= deg2rad;

	var  e  = 0.016708634 + T * ( -0.000042037 - 0.0000001267 * T );

	var  y  = calOblique( JD );
		 y  = tan4deg( y / 2.0 );
		 y  = y * y;

	var  E  = y * Math.sin( 2 * L0 ) - 2.0 * e * Math.sin( M );
		 E += 4.0 * e * y * Math.sin( M ) * Math.cos( 2.0 * L0 );
		 E -= y * y * Math.sin( 4.0 * L0 ) / 2.0;
		 E -= 5.0 * e * e * Math.sin( 2.0 * M ) / 4.0;

	E /= deg2rad;
	return ( E * 4.0 );
}

////////////////////////////////

// カレンダー関係。
function calDayOfWeek(year, month, day){
	var JD = calJDz(year, month, day);
	var you = (JD + 1) % 7;
	return you;
}

function chkLeap( year ){
	var chk = 0;

	if(year %   4 == 0) chk = 1;
	if(year % 100 == 0) chk = 0;
	if(year % 400 == 0) chk = 1;

	return chk;
}

function maxday(year, month){
	var mday = [31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
	var md = mday[ month ];
	if( month == 2 && chkLeap( year ) ){
		md = 29;
	}

	return md;
}

// ΔＴを管理する関数
// formula A : Notes Scientifiques et Techniques du Bureau des Longitudes, nr. S055
// from ftp://cyrano-se.obspm.fr/pub/6_documents/4_lunar_tables/newexp.pdf
// formula B : Polynomial Expressions for Delta T (ΔT)
// from https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
// formula C : Delta T : Polynomial Approximation of Time Period 1620-2013
// from https://www.hindawi.com/archive/2014/480964/ (license: CC-BY-3.0)
function correctTDT(JD){
	var year = ( JD - 2451545.0 ) / 365.2425 + 2000.0;
	var t, dt;

	if( year < 948 ){ // formula A.26
		t = ( JD - 2451545.0 ) / 36525.0;
		dt = 2177.0 + t * ( 497.0 + t * 44.1 );
	} else if( year < 1600 ){ // formula A.25
		t = ( JD - 2451545.0 ) / 36525.0;
		dt =  102.0 + t * ( 102.0 + t * 25.3 );
	} else if( year < 1620 ){ // formula B
		t = year - 1600;
		dt = 120 + t * ( -0.9808 + t * ( -0.01532 + t / 7129 ) );
	} else if( year < 2014 ){ // formula C
		// last elements are sentinels.
		var tep = [     1620,     1673,     1730,      1798,      1844,     1878,      1905,      1946,      1990,  2014 ];
		var tk  = [    3.670,    3.120,    2.495,     1.925,     1.525,    1.220,     0.880,     0.455,     0.115, 0.000 ];
		var ta0 = [   76.541,   10.872,   13.480,    12.584,     6.364,   -5.058,    13.392,    30.782,    55.281, 0.000 ];
		var ta1 = [ -253.532,  -40.744,   13.075,     1.929,    11.004,   -1.701,   128.592,    34.348,    91.248, 0.000 ];
		var ta2 = [  695.901,  236.890,    8.635,    60.896,   407.776,  -46.403,  -279.165,    46.452,    87.202, 0.000 ];
		var ta3 = [-1256.982, -351.537,   -3.307, -1432.216, -4168.394, -866.171, -1282.050,  1295.550, -3092.565, 0.000 ];
		var ta4 = [  627.152,   36.612, -128.294,  3129.071,  7561.686, 5917.585,  4039.490, -3210.913,  8255.422, 0.000 ];

		var i = 0;
		for( var j = 0; j < tep.length; j++ ){
			if( tep[ j ] <= year && year < tep[ j + 1 ] ){
				i = j;
				break;
			}
		}
		var k  = tk[ i ];
		var a0 = ta0[ i ];
		var a1 = ta1[ i ];
		var a2 = ta2[ i ];
		var a3 = ta3[ i ];
		var a4 = ta4[ i ];

		var u = k + ( year - 2000 ) / 100;
		dt = a0 + u * ( a1 + u * ( a2 + u * ( a3 + u * a4 ) ) );
	} else { // formula A.25
		t = ( JD - 2451545.0 ) / 36525.0;
		dt =  102.0 + t * ( 102.0 + t * 25.3 );
		if( year < 2100 ){
			dt += 0.37 * ( year - 2100 ); // from "Astronomical Algorithms" p.78
		}
	}

	dt /= 86400.0;
	return dt;
}

function advanceDate( date, step ){
	return encodeDate(cnvCalendar(calJDz(decodeDate(date)) + step));
}

function calDist(sy, sm, sd, ey, em, ed){
	return calJDz(ey, em, ed) - calJDz(sy, sm, sd);
}

function decodeDate( date ){
	var ye = Math.floor( date / 10000 );
	var mo = Math.floor((date % 10000) / 100);
	var da =             date %   100;

	var res = new Array(ye, mo, da);
	return res;
}

function decodeTime( time ){
	var ho = Math.floor(time / 100);
	var mi = fmod(time, 100);

	var res = new Array(ho, mi);
	return res;
}

function encodeDate(ye, mo, da){
	return ye * 10000 + mo * 100 + da;
}

function encodeTime(ho, mi){
	return ho * 100 + mi;
}
