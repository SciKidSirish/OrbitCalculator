import numpy as np
import astropy
import datetime
from astropy.coordinates import  Galactic, FK4, FK5, AltAz  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame, EarthLocation, AltAz, ICRS
from astropy.time import Time
import math



tList = []

parallaxDone = False
for i in range(1,11):
  print("Observation {}:".format(i))
  t = input("What is the date(in YYYY-MM-DD format with no spaces)?  ") + "T" + input("What is the time in HH:MM:SS 24hr format?   ") + "Z"
  alt = float(input("What is the altitude in degrees (exlude symbol)?   "))
  az = float(input("What is the azimuth in degrees (exlude symbol)?   "))
  lat = float(input("What is the lattitude in degrees (exlude symbol)?   "))
  lon = float(input("What is the longitude in degrees (exlude symbol)?   "))
  if parallaxDone == False:
    if input("Would you like to add a parallax measurement (y/n)?   ") == "y":
      print("Parallax measurements:")
      pt = input("What is the date(in YYYY-MM-DD format with no spaces)?  ") + "T" + input("What is the time in HH:MM:SS 24hr format?   ") + "Z"
      palt = float(input("What is the altitude in degrees (exlude symbol)?   "))
      paz = float(input("What is the azimuth in degrees (exlude symbol)?   "))
      plat = float(input("What is the lattitude in degrees (exlude symbol)?   "))
      plon = float(input("What is the longitude in degrees (exlude symbol)?   "))
      p_angle = (180-(paz+az))/2
      cLength = (2* 0.0000426354 * math.sin(math.radians(abs(plon-lon))))/2
      eDist = math.asin(math.radians(p_angle) % 1) * cLength
      gcrsCoord = astropy.coordinates.GCRS(ra = paz * u.deg, dec = (palt+90-paz) % 90 * u.deg)
      point = gcrsCoord.transform_to(ICRS)
      point = str(point)
      pre = "<ICRS Coordinate: (ra, dec) in deg"
      point = point.removeprefix(pre)
      point = point.lstrip()
      ind = point.index('>')
      point = point[:ind]
      point = point.strip('()')
      point = [cLength, list(map(float, point.split(', ')))[0]]
      parallaxDone = True



  # path = input("What is the path to your image?")
  # bList.append(brightness(path))
  obstime = Time(t)
  location = EarthLocation.from_geodetic(18 * u.deg, lat=0 * u.deg, height=100 * u.m)
  altaz_frame = AltAz(obstime=obstime, location=location)
  target = SkyCoord(alt=alt * u.deg, az=az * u.deg, frame=altaz_frame)
  dirs_altaz_offset = SkyCoord(
      lon=[-0.02, 0.01, 0.0, 0.0, 0.0] * u.rad,
      lat=[0.0, 0.2, 0.0, -0.3, 0.1] * u.rad,
      frame=target.skyoffset_frame()
  )
  dirs_altaz = dirs_altaz_offset.transform_to(altaz_frame)
  dirs_icrs = dirs_altaz.transform_to(ICRS())
  target_icrs = target.transform_to(ICRS())
  dirs_icrs_offset = dirs_icrs.transform_to(target_icrs.skyoffset_frame())
  eq = str(dirs_icrs_offset)
  pre = "<SkyCoord (SkyOffsetICRS: rotation=0.0 deg, origin=<ICRS Coordinate: (ra, dec) in deg"
  eq = eq.removeprefix(pre)
  eq = eq.lstrip()
  ind = eq.index('>')
  eq = eq[:ind]
  eq = eq.strip('()')
  eqList = list(map(float, eq.split(', ')))

  date_format = datetime.datetime.strptime(t,
                                          "%Y-%m-%dT%H:%M:%SZ")
  unix_time = datetime.datetime.timestamp(date_format)
  eqList.append(unix_time/31536000)
  tList.append(eqList)

  
  
  import astropy
import astropy.coordinates as c
import astropy.units as u
from statistics import mean
from sympy import Plane, Point3D
import numpy as np
import math

radec = tList

def cstrip(cart):
  cart = str(cart)
  cart = cart.replace("<Quantity ", "")
  cart = cart.replace(">", "")
  cart = cart.replace("(", "")
  cart = cart.replace(")", "")
  cart = cart.split(", ")
  for i in range(0, len(cart)):
    cart[i] = float(cart[i])
  return cart

angList = []

for i in range(0,8):
   cart1 = cstrip(c.spherical_to_cartesian(1, radec[i][0], radec[i][1]))
   cart2 = cstrip(c.spherical_to_cartesian(1, radec[i+1][0], radec[i+1][1]))
   cart3 = cstrip(c.spherical_to_cartesian(1, radec[i+2][0], radec[i+2][1]))
   A = Plane(Point3D(cart1[0], cart1[1], cart1[2]), Point3D(cart2[0], cart2[1], cart2[2]), Point3D(cart3[0], cart3[1], cart3[2]))
   E = Plane(Point3D(1, 0, 0), Point3D(-1, 0, 0), Point3D(0, 1, 0))
  #  print(A.angle_between(E))
   angList.append(float(A.angle_between(E)))

# print(angList)

avgAng = mean(angList)


#----------------------------------------------------Ellipse
ecoord=[]

for i in range(0, len(radec)):
  ecoord.append([radec[i][0], radec[i][2]])

P = (abs(ecoord[0][1]-ecoord[1][1])/abs(ecoord[0][0]-ecoord[1][0]))*(360 * u.deg)

sma = np.cbrt(P**2)

# print(sma)

E = ((-1*(point[0]*math.cos(point[1])))+(math.sqrt( (point[0]*math.cos(point[1])) ** 2 + 4*(sma**2) - 4*sma*point[0])))/(2*sma)


print("Here are the characteristics of the orbit: Semi-major Axis: {} | Eccentricity: {} | Angle from plane of equator: {}".format(sma, E, avgAng))
