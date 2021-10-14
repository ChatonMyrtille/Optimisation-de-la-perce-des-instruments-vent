import math
import cmath
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import random as rnd
import copy

# density of air
rho = 1.225

# viscosity of air
eta = 0.0000181
    
# cylinder cell
class CylCell:

    # create a cylinder
    def __init__(self, a, L):
    
        # radius of the cylinder
        self.a = a
        
        # length of the cylinder
        self.L = L
    
    # get the transfer matrix
    def transferMatrix(self, f):
        return Tcyl(self.a, self.L, f)

# tone hole cell
class HoleCell:

    # create a tone hole
    def __init__(self, a, b, h):
    
        # radius of the cylinder
        self.a = a
        
        # radius of the hole
        self.b = b
        
        # length of the chimney
        self.h = h
    
    # get the transfer matrix
    def transferMatrix(self, f, open):
        return Thole(self.a, self.b, self.h, f, open)
    
# cone cell
class ConeCell:

    # create a cone cell
    def __init__(self, a1, a2, L):
    
        # input radius
        self.a1 = a1
        
        # output radius
        self.a2 = a2
        
        # lenght between input plane and output plane
        self.L = L
    
    # get the transfer matrix
    def transferMatrix(self, f):
        return Tcone(self.a1, self.a2, self.L, f)

# instrument
class Instrument:

    def __init__(self):
        self.cells = []
        self.impedanceResonances = {}
    
    def addCell(self, cell):
        self.cells.append(cell)
    
    # get the impedance
    def impedance(self, f, fingering):
    
        # initialize current tone hole
        currentToneHole = 1
        
        # transfer matrix
        T = np.identity(2, dtype=complex)
        
        # for each cell of the instrument
        for cell in self.cells:
                
            # when the cell is a tone hole
            if cell.__class__.__name__ == "HoleCell":
                
                # when the hole is open
                if fingering[currentToneHole-1] == 1:
                    
                    # update the transfer matrix
                    T = np.dot(T, cell.transferMatrix(f, 1))
                        
                # when the hole is closed
                else:
                    
                    # update the transfer matrix
                    T = np.dot(T, cell.transferMatrix(f, 0))
                        
                # increment the current tone hole
                currentToneHole += 1
                
            # when the cell is a cylinder or a cone
            else:
            
                #update the transfer matrix
                T = np.dot(T, cell.transferMatrix(f))
        
        # angular frequency
        omega = 2*constants.pi*f
        
        # wave number
        k = omega/343.0
        
        # output radius
        a = self.cells[-1].a2
        
        # cross section area
        S = constants.pi*a**2
                
        # radiation impedance
        Zr = ((rho*343.0)/S) * complex(0.25*(k*a)**2, 0.6*k*a)
        T = np.dot(T,np.matrix([[Zr],[complex(1,0)]], dtype=complex))
        
        # input impedance
        Zin = T.item(0)/T.item(1)
        return Zin
    
    # draw the instrument
    def draw(self):
    
        # create a new image
        image = Image.new("RGBA", (1000, 500), (255, 255, 255, 255))
        draw = ImageDraw.Draw(image)
        fontSmall = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 10)
        fontLarge = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 15)
        white = (255,255,255,255)
        black = (0,0,0,255)
        
        # draw cones
        currentLength = 100 - self.cells[0].L*1000
        for i in range(len(self.cells)):

            # draw cones
            if self.cells[i].__class__.__name__ == "ConeCell":
            
                # draw first cone
                if i == 0:
                    
                    draw.line((currentLength + self.cells[i].L*1000, 250 - 0.0045*1000,  + currentLength + self.cells[i].L*1000, 250 - 0.0095*1000), fill=black)
                    
                    draw.line((currentLength + self.cells[i].L*1000, 250 + 0.0045*1000,  + currentLength + self.cells[i].L*1000, 250 + 0.0095*1000), fill=black)
                    
                    draw.line((currentLength, 250 - self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 - self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
            
                # draw second cone
                elif i == 1:
                    
                    draw.line((currentLength, 250 - self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 - self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength, 250 - (self.cells[i].a1+0.0075)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i].a2+0.0075)*1000), fill=black)
                    
                    draw.line((currentLength, 250 + (self.cells[i].a1+0.0075)*1000, currentLength + self.cells[i].L*1000, 250 + (self.cells[i].a2+0.0075)*1000), fill=black)
                    
                # draw third cone
                elif i == 2:
                
                    draw.line((currentLength, 250 - self.cells[i].a1*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength, 250 - (self.cells[i].a1+0.0075)*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - (self.cells[i].a2+self.cells[i+1].h)*1000), fill=black)
                    
                    draw.line((currentLength, 250 + (self.cells[i].a1+0.0075)*1000, currentLength + (self.cells[i].L)*1000, 250 + (self.cells[i].a2+0.0075)*1000), fill=black)
                
                # draw last cone
                elif i == len(self.cells)-1:
                
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - self.cells[i].a1*1000, currentLength + (self.cells[i].L)*1000, 250 - self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                    
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - (self.cells[i].a1+self.cells[i-1].h)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i].a2+0.0035)*1000), fill=black)
                    
                    draw.line((currentLength, 250 + (self.cells[i].a1+0.0075)*1000, currentLength + (self.cells[i].L)*1000, 250 + (self.cells[i].a2+0.0075)*1000), fill=black)
                    
                # draw other cones
                else:
                    
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - self.cells[i].a1*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - self.cells[i].a2*1000), fill=black)
                
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - (self.cells[i].a1+self.cells[i-1].h)*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - (self.cells[i].a2+self.cells[i+1].h)*1000), fill=black)
                
                    draw.line((currentLength, 250 + (self.cells[i].a1+0.0075)*1000, currentLength + (self.cells[i].L)*1000, 250 + (self.cells[i].a2+0.0075)*1000), fill=black)
                    
                # increment current length
                currentLength += self.cells[i].L*1000
        
        # draw holes
        currentLength = 100 - self.cells[0].L*1000
        for cell in self.cells:
            
            if cell.__class__.__name__ == "ConeCell":
                currentLength += cell.L*1000
            elif cell.__class__.__name__ == "HoleCell":
                
                # draw chimney
                draw.line((currentLength - cell.b*1000, 250 - (cell.a+cell.h)*1000, currentLength - cell.b*1000, 250 - cell.a*1000), fill=black)
                draw.line((currentLength + cell.b*1000, 250 - (cell.a+cell.h)*1000, currentLength + cell.b*1000, 250 - cell.a*1000), fill=black)
            
        # draw length axis
        for i in range(0, 800, 100):
            w, h = draw.textsize(str(i)+" mm", font=fontLarge)
            if i <= 600:
                draw.line((100 + i, 300, 100 + i + 100, 300), fill=black)
                draw.line((100 + i + 25, 305, 100 + i + 25, 295), fill=black)
                draw.line((100 + i + 50, 305, 100 + i + 50, 295), fill=black)
                draw.line((100 + i + 75, 305, 100 + i + 75, 295), fill=black)
            draw.line((100 + i, 305, 100 + i, 295), fill=black)
            draw.text((100 + i - w/2, 315 - h/2), str(i)+" mm", fill=black, font=fontLarge)
        
        image.show()
    
    # plot reactance
    def plotReactance(self, fingering):
        # plotting data
        xdata = []
        ydata = []

        for f in range(100, 2500):

            # impedance of the instrument
            Zin = self.impedance(f, fingering)
            
            # add data
            xdata.append(f)
            ydata.append(Zin.imag)

        # plot the reactance of the instrument
        plt.figure(figsize=(10,5))
        plt.plot(xdata,ydata,"black")
        plt.title("Réactance de l'instrument")
        plt.xlabel('fréquence (Hz)')
        plt.ylabel('réactance (Pa·s.m-3)')
        plt.show()
    
    # plot impedance
    def plotImpedance(self, fingering):
        # plotting data
        xdata = []
        ydata = []

        for f in range(100, 2500):

            # impedance of the instrument
            Zin = self.impedance(f, fingering)
            
            # add data
            xdata.append(f)
            ydata.append(abs(Zin))

        # plot the impedance of the instrument
        plt.figure(figsize=(10,5))
        plt.yscale("log")
        plt.plot(xdata,ydata,"black")
        plt.title("Impédance de l'instrument")
        plt.xlabel('fréquence (Hz)')
        plt.ylabel('impédance (Pa·s.m-3)')
        plt.show()
    
    # get the playing frequency of the instrument
    def playingFrequencyMaxima(self, fingering, resonances):
    
        approxFrequency = 0.0
        flag = True
        
        for f in range(100, 3000):
            Zin = self.impedance(f, fingering)
                        
            # sign is negative
            if np.sign(Zin.imag) == -1:
                if flag:
                    if resonances == 1:
                        approxFrequency = f
                        break
                    else:
                        resonances = resonances - 1
                        flag = False
            else:
                flag = True
                                
        for f in np.arange(approxFrequency-1,approxFrequency+1,0.01):

            # impedance of the instrument
            Zin = self.impedance(f, fingering)
            #print("f=" + "%.2f" % round(f,2) + ", Zin(ω)=" + str(Zin.imag))
            
            # return playing frequency
            if np.sign(Zin.imag) == -1:
                return f
                
    # get the playing frequency of the instrument
    def playingFrequencyMinima(self, fingering, resonances):
    
        approxFrequency = 0.0
        flag = False
        
        for f in range(100, 5000):
            Zin = self.impedance(f, fingering)
                        
            # sign is positive
            if np.sign(Zin.imag) == +1:
                if flag:
                    if resonances == 1:
                        approxFrequency = f
                        break
                    else:
                        resonances = resonances - 1
                        flag = False
            else:
                flag = True
                                
        for f in np.arange(approxFrequency-1,approxFrequency+1,0.01):

            # impedance of the instrument
            Zin = self.impedance(f, fingering)
            #print("f=" + "%.2f" % round(f,2) + ", Zin(ω)=" + str(Zin.imag))
            
            # return playing frequency
            if np.sign(Zin.imag) == +1:
                return f
                
    # evolve the instrument
    def evolve(self, value):
    
        # change input diameter
        #d = rnd.random()*0.001 - 0.0005
        #self.cells[1].a1 += d
        #if self.cells[1].a1 < 0.0050: self.cells[1].a1 = 0.0050
        #if self.cells[1].a1 > 0.0100: self.cells[1].a1 = 0.0100
        #self.cells[1].a2 += d
        #if self.cells[1].a2 < 0.0050: self.cells[1].a2 = 0.0050
        #if self.cells[1].a2 > 0.0100: self.cells[1].a2 = 0.0100
        #self.cells[2].a1 = self.cells[1].a2
        
        # change input length
        self.cells[1].L += rnd.random()*value - value/2.0
        if self.cells[1].L < 0.125: self.cells[1].L = 0.125
        if self.cells[1].L > 0.175: self.cells[1].L = 0.175
        
        for i in range(len(self.cells)):
        
            if self.cells[i].__class__.__name__ == "HoleCell":
            
                # change bore radius
                self.cells[i].a += rnd.random()*value - value/2.0
                if self.cells[i].a < 0.005: self.cells[i].a = 0.005
                if self.cells[i].a > 0.010: self.cells[i].a = 0.010
                self.cells[i-1].a2 = self.cells[i].a
                self.cells[i+1].a1 = self.cells[i].a
                
                # change hole radius
                self.cells[i].b += rnd.random()*value - value/2.0
                
                if self.cells[i].b < 0.001: self.cells[i].b = 0.001
                if self.cells[i].b > 0.0055: self.cells[i].b = 0.0055
                
                # change chimney length
                self.cells[i].h += rnd.random()*value - value/2.0
                if self.cells[i].h < 0.001: self.cells[i].h = 0.001
                if self.cells[i].h > 0.005: self.cells[i].h = 0.005
                
        for i in range(2, len(self.cells)):
            
            if self.cells[i].__class__.__name__ == "ConeCell":
                
                # change hole positions
                self.cells[i].L += rnd.random()*value - value/2.0
                if self.cells[i].L < 0.020: self.cells[i].L = 0.020
                if self.cells[i].L > 0.100: self.cells[i].L = 0.100
            
        # change output diameter
        self.cells[-1].a2 += rnd.random()*value - value/2.0
        if self.cells[-1].a2 < 0.005: self.cells[-1].a2 = 0.005
        if self.cells[-1].a2 > 0.010: self.cells[-1].a2 = 0.010
            

    # plot distance in Hz
    def plotDistanceHz(self):
        x = ["Ré3", "Ré3#", "Mi3", "Fa3", "Fa3#", "Sol3", "Sol3#", "La3", "La3#", "Si3", "Do4", "Do4#", "Ré4", "Ré4#", "Mi4", "Fa4", "Fa4#", "Sol4", "Sol4#", "La4", "La4#", "Si4", "Do5", "Do5#", "Ré5", "Ré5#", "Mi5", "Fa5", "Fa5#", "Sol5", "Sol5#", "La5"]
        y1 = []
        y2 = []
        y3 = []

        for note in fluteNotes:
            tf = theoreticalFrequency[note]
            f1 = self.impedanceResonances[note][0]
            f2 = self.impedanceResonances[note][1]
        
            # x.append(note)
            y1.append(tf)
            y2.append(f1)
            y3.append(f2)
                    
        plt.figure(figsize=(15,10))
        plt.plot(x, y1, "-", color="grey", label="fréquence théorique")
        plt.plot(x, y2, "+--", color="black", label="minima de l'impédance")
        plt.plot(x, y3, "+--", color="black")
        plt.legend(loc="upper left")
        plt.title("Distance entre la fréquence de l'instrument et la fréquence théorique")
        plt.xlabel("notes")
        plt.ylabel("fréquence (Hz)")
        plt.show()
    
    # plot distance in cents
    def plotDistanceCents(self):
        x = ["Ré3", "Ré3#", "Mi3", "Fa3", "Fa3#", "Sol3", "Sol3#", "La3", "La3#", "Si3", "Do4", "Do4#", "Ré4", "Ré4#", "Mi4", "Fa4", "Fa4#", "Sol4", "Sol4#", "La4", "La4#", "Si4", "Do5", "Do5#", "Ré5", "Ré5#", "Mi5", "Fa5", "Fa5#", "Sol5", "Sol5#", "La5"]
        y = []
        z = []

        for note in fluteNotes:
            # x.append(note)
            z.append(0)
            
            f1 = self.impedanceResonances[note][0]
            f2 = self.impedanceResonances[note][1]
            
            tf = theoreticalFrequency[note]
            
            if note in fluteFirstRegister:
                a = tf
                b = f1
                y.append(1200*np.log2(b/a))
                
            elif note in fluteSecondRegister:
                a = tf
                b = f2
                y.append(1200*np.log2(b/a))
                    
            elif note in fluteThirdRegister:
                y.append(float("nan"))

        plt.figure(figsize=(15,10))
        plt.plot(x, y, "+--", color="black", label="fréquence de l'instrument")
        plt.plot(x, z, "-", color="grey", label="fréquence théorique")
        plt.legend(loc="upper left")
        plt.title("Distance entre la fréquence de l'instrument et la fréquence théorique")
        plt.xlabel("notes")
        plt.ylabel("distance (cents)")
        plt.show()
    
    # print the parameters
    def printParameters(self):
        currentLength = -self.cells[0].L
        holeNumber = 1
        print("mouth piece diameter: ø " + "%.2f" % round(self.cells[0].a1*2*1000,2) + " mm")
        print("mouth piece length: " + "%.2f" % round(self.cells[0].L*1000,2) + " mm")
        print("first cylinder diameter: ø " + "%.2f" % round(self.cells[1].a1*2*1000,2) + " mm")
        print("first cylinder length: " + "%.2f" % round(self.cells[1].L*1000,2) + " mm")
        print("hole # \t\tposition \tbore diameter \thole diameter \tchimney length")
        for i in range(len(self.cells)):
            if self.cells[i].__class__.__name__ == "CylCell" or self.cells[i].__class__.__name__ == "ConeCell":
                currentLength += self.cells[i].L
            if self.cells[i].__class__.__name__ == "HoleCell":
                if holeNumber == 1:
                    print(str(holeNumber) + "st hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].a*2*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                elif holeNumber == 2:
                    print(str(holeNumber) + "nd hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].a*2*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                elif holeNumber == 3:
                    print(str(holeNumber) + "rd hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].a*2*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                else:
                    print(str(holeNumber) + "th hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].a*2*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                holeNumber += 1
        print("output diameter: ø " + "%.2f" % round(self.cells[-1].a2*2*1000,2) + " mm")
        print("last cone length: " + "%.2f" % round(self.cells[-1].L*1000,2) + " mm")

# transfer matrix for a cylinder
def Tcyl(a, L, f):

    # angular frequency
    omega = 2*constants.pi*f

    # wave number
    k = omega/343.0
    
    rv = a * math.sqrt(rho*omega/eta)
    Z0 = (rho * 343.0) / (constants.pi * a**2)
    Zc = Z0 * complex(1 + 0.369/rv, -(0.369/rv + 1.149/rv**2))
    gamma = k * complex(1.045/rv + 1.080/rv**2, 1 + 1.045/rv)
    return np.matrix([[cmath.cosh(gamma*L),Zc*cmath.sinh(gamma*L)],[(1.0/Zc)*cmath.sinh(gamma*L),cmath.cosh(gamma*L)]], dtype=complex)

# transfer matrix for a tone hole
def Thole(a, b, h, f, open):

    # angular frequency
    omega = 2*constants.pi*f

    # wave number
    k = omega/343.0

    delta = b/a
    ta = b * (-0.37 + 0.078*delta) * delta**2
    ts = b * (0.82 - 0.193*delta - 1.09*delta**2 + 1.27*delta**3 - 0.71*delta**4)
    ma = (rho*ta)/(constants.pi*a**2)
    ms = (rho*ts)/(constants.pi*b**2)
    
    # series impedance
    Za = complex(0, omega*ma)
    
    # shunt impedance
    Zs = complex(0, omega*ms)
    
    # cross section area
    Sh = constants.pi*b**2
    
    # radiation impedance
    ZL = ((rho*343.0)/(Sh))*complex(0.5*(k*a)**2, 0.82*k*a)
        
    # input impedance
    Zh = ((rho*343.0)/Sh) * (ZL + complex(0,((rho*343.0)/Sh)*math.tan(k*h)))/((rho*343.0)/Sh + ZL*complex(0,math.tan(k*h)))
    
    if open == 0:
        Zh = (rho*343.0**2)/(omega*Sh*h)
    
    # total shunt impedance
    Zst = Zs - Za/4 + Zh
    
    return np.matrix([[1.0,Za/2.0],[0.0,1.0]], dtype=complex).dot(np.matrix([[1.0,0.0],[1.0/Zst,1.0]], dtype=complex)).dot(np.matrix([[1.0,Za/2],[0.0,1.0]], dtype=complex))

# transfer matrix for a cone
def Tcone(a1, a2, L, f):

    # input diameter and output diameter are equals
    if a1 == a2:
    
        # consider the cone like a cylinder
        return Tcyl(a1, L, f)

    # angular frequency
    omega = 2*constants.pi*f

    # wave number
    k = omega/343.0

    ax = -L/(-a2+a1)
    ay = -L/(+a2-a1)
    bx = +a2*ax
    by = -a2*ay
    x = -(by-bx)/(ay-ax)
    y = ax*x+bx
    x1 = math.sqrt((y-L)**2)
    x2 = math.sqrt(y**2)
    Zc = (rho*343.0)/(constants.pi*a1*a2)
    aeq = L*(a1/x1)*(1.0/np.log(1+L/x1))
    gamma = complex(aeq, omega/343.0 + aeq)
    kc = complex(0, -gamma)
    return np.matrix([[(a2/a1)*cmath.cos(kc*L)-cmath.sin(kc*L)/(k*x1),complex(0,Zc*cmath.sin(kc*L))],[(1.0/Zc)*(complex(0, (1+(1.0/((k**2)*x1*x2)))*cmath.sin(kc*L))+complex(0,-(1.0/x1-1.0/x2)*cmath.cos(kc*L)/k)),(a1/a2)*cmath.cos(kc*L)+cmath.sin(kc*L)/(k*x2)]], dtype=complex)

# creat an baroque flute
flute = Instrument()

# mouthpiece (16.5 mm)
flute.addCell(ConeCell(0.0045, 0.0045, 0.0165))

# bore (150 mm)
flute.addCell(ConeCell(0.0095, 0.0095, 0.150))

# bore (80 mm)
flute.addCell(ConeCell(0.0095, 0.0090, 0.080))

# 1st tone hole
flute.addCell(HoleCell(0.0090, 0.0030, 0.0035))

# bore (35 mm)
flute.addCell(ConeCell(0.0090, 0.00875, 0.035))

# 2nd tone hole
flute.addCell(HoleCell(0.00875, 0.0030, 0.0035))

# bore (40 mm)
flute.addCell(ConeCell(0.00875, 0.0085, 0.040))

# 3rd tone hole
flute.addCell(HoleCell(0.0085, 0.0025, 0.0035))

# bore (60 mm)
flute.addCell(ConeCell(0.0085, 0.0080, 0.060))

# 4th tone hole
flute.addCell(HoleCell(0.0080, 0.0030, 0.0035))

# bore (35 mm)
flute.addCell(ConeCell(0.0080, 0.0075, 0.035))

# 5th tone hole
flute.addCell(HoleCell(0.0075, 0.0030, 0.0035))

# bore (40 mm)
flute.addCell(ConeCell(0.0075, 0.00725, 0.040))

# 6th tone hole
flute.addCell(HoleCell(0.00725, 0.0020, 0.0035))

# bore (60 mm)
flute.addCell(ConeCell(0.00725, 0.00625, 0.060))

# 7th tone hole
flute.addCell(HoleCell(0.00625, 0.0030, 0.0035))

# bore (55 mm)
flute.addCell(ConeCell(0.00625, 0.0070, 0.055))

# fingering chart for a 7 tone holes baroque flute
fingeringChartFlute = {}
fingeringChartFlute["D4"] = [0,0,0,0,0,0,0]
fingeringChartFlute["D#4"] = [0,0,0,0,0,0,1]
fingeringChartFlute["E4"] = [0,0,0,0,0,1,0]
fingeringChartFlute["F4"] = [0,0,0,0,1,0,0]
fingeringChartFlute["F#4"] = [0,0,0,0,1,1,1]
fingeringChartFlute["G4"] = [0,0,0,1,1,1,1]
fingeringChartFlute["G#4"] = [0,0,1,0,0,0,1]
fingeringChartFlute["A4"] = [0,0,1,1,1,1,1]
fingeringChartFlute["A#4"] = [0,1,0,0,0,1,0]
fingeringChartFlute["B4"] = [0,1,1,1,1,1,1]
fingeringChartFlute["C5"] = [1,0,0,1,1,1,1]
fingeringChartFlute["C#5"] = [1,1,1,1,1,1,0]
fingeringChartFlute["D5"] = [1,0,0,0,0,0,0]
fingeringChartFlute["D#5"] = [1,0,0,0,0,0,1]
fingeringChartFlute["E5"] = [0,0,0,0,0,1,0]
fingeringChartFlute["F5"] = [0,0,0,0,1,0,0]
fingeringChartFlute["F#5"] = [0,0,0,0,1,1,1]
fingeringChartFlute["G5"] = [0,0,0,1,1,1,1]
fingeringChartFlute["G#5"] = [0,0,1,0,1,1,1]
fingeringChartFlute["A5"] = [0,0,1,1,1,1,1]
fingeringChartFlute["A#5"] = [0,1,0,1,1,1,1]
fingeringChartFlute["B5"] = [0,1,1,1,1,1,1]
fingeringChartFlute["C6"] = [1,0,1,0,0,0,1]
fingeringChartFlute["C#6"] = [1,0,0,0,1,1,1]
fingeringChartFlute["D6"] = [1,0,0,1,1,1,1]
fingeringChartFlute["D#6"] = [0,0,0,1,0,0,1]
fingeringChartFlute["E6"] = [0,0,1,1,0,0,1]
fingeringChartFlute["F6"] = [0,0,1,0,1,1,1]
fingeringChartFlute["F#6"] = [0,0,1,0,1,1,0]
fingeringChartFlute["G6"] = [0,1,0,1,1,1,0]
fingeringChartFlute["G#6"] = [1,1,0,1,1,0,1]
fingeringChartFlute["A6"] = [1,0,0,0,0,1,0]

A440 = 440
A415 = 415
r = A415/A440

# theoretical frequencies
theoreticalFrequency = {}
theoreticalFrequency["D3"] = 146.8324 * r
theoreticalFrequency["D#3"] = 155.5635 * r
theoreticalFrequency["E3"] = 164.8138 * r
theoreticalFrequency["F3"] = 174.6141 * r
theoreticalFrequency["F#3"] = 184.9972 * r
theoreticalFrequency["G3"] = 195.9977 * r
theoreticalFrequency["G#3"] = 207.6523 * r
theoreticalFrequency["A3"] = 220.0000 * r
theoreticalFrequency["A#3"] = 233.0819 * r
theoreticalFrequency["B3"] = 246.9417 * r
theoreticalFrequency["C4"] = 261.6256 * r
theoreticalFrequency["C#4"] = 277.1826 * r
theoreticalFrequency["D4"] = 293.6648 * r
theoreticalFrequency["D#4"] = 311.1270 * r
theoreticalFrequency["E4"] = 329.6276 * r
theoreticalFrequency["F4"] = 349.2282 * r
theoreticalFrequency["F#4"] = 369.9944 * r
theoreticalFrequency["G4"] = 391.9954 * r
theoreticalFrequency["G#4"] = 415.3047 * r
theoreticalFrequency["A4"] = 440.0000 * r
theoreticalFrequency["A#4"] = 466.1638 * r
theoreticalFrequency["B4"] = 493.8833 * r
theoreticalFrequency["C5"] = 523.2511 * r
theoreticalFrequency["C#5"] = 554.3653 * r
theoreticalFrequency["D5"] = 587.3295 * r
theoreticalFrequency["D#5"] = 622.2540 * r
theoreticalFrequency["E5"] = 659.2551 * r
theoreticalFrequency["F5"] = 698.4565 * r
theoreticalFrequency["F#5"] = 739.9888 * r
theoreticalFrequency["G5"] = 783.9909 * r
theoreticalFrequency["G#5"] = 830.6094 * r
theoreticalFrequency["A5"] = 880.0000 * r
theoreticalFrequency["A#5"] = 932.3275 * r
theoreticalFrequency["B5"] = 987.7666 * r
theoreticalFrequency["C6"] = 1046.502 * r
theoreticalFrequency["C#6"] = 1108.731 * r
theoreticalFrequency["D6"] = 1174.659 * r
theoreticalFrequency["D#6"] = 1244.508 * r
theoreticalFrequency["E6"] = 1318.510 * r
theoreticalFrequency["F6"] = 1396.913 * r
theoreticalFrequency["F#6"] = 1479.978 * r
theoreticalFrequency["G6"] = 1567.982 * r
theoreticalFrequency["G#6"] = 1661.219 * r
theoreticalFrequency["A6"] = 1760.000 * r

# range of the baroque flute
fluteNotes = ["D4", "D#4", "E4", "F4", "F#4", "G4", "G#4", "A4", "A#4", "B4", "C5", "C#5", "D5", "D#5", "E5", "F5", "F#5", "G5", "G#5", "A5", "A#5", "B5", "C6", "C#6", "D6", "D#6", "E6", "F6", "F#6", "G6", "G#6", "A6"]

# flute first register
fluteFirstRegister = ["D4", "D#4", "E4", "F4", "F#4", "G4", "G#4", "A4", "A#4", "B4", "C5", "C#5"]

# flute second register
fluteSecondRegister = ["D5", "D#5", "E5", "F5", "F#5", "G5", "G#5", "A5", "A#5", "B5"]

# flute third register
fluteThirdRegister = ["C6", "C#6", "D6", "D#6", "E6", "F6", "F#6", "G6", "G#6", "A6"]

# the distance from theoretical frequency for each notes
distances = []

# for each note
for note in fluteNotes:

    # get the playing frequency
    f1 = flute.playingFrequencyMinima(fingeringChartFlute[note], 1)
    f2 = flute.playingFrequencyMinima(fingeringChartFlute[note], 2)
    #print(note + ": " + "%.2f" % round(f1,2) + " Hz, " + "%.2f" % round(f2,2) + " Hz")
    
    flute.impedanceResonances[note] = (f1,f2)
    
    # theoretical frequency
    tf = theoreticalFrequency[note]
    
    if note in fluteFirstRegister:
        distances.append(abs(1200*np.log2(f1/tf)))

    elif note in fluteSecondRegister:
        distances.append(abs(1200*np.log2(f2/tf)))

    elif note in fluteThirdRegister:
        distances.append(0)
        
# total distance
totalDistance = sum(distances)

# worst distance
worstDistance = max(distances)
    
print("total distance from theoretical frequency: " + "%.2f" % round(totalDistance,2) + " cents")
print("worst distance from theoretical frequency: " + "%.2f" % round(worstDistance,2) + " cents")

# number of iterations
number_of_iterations = 0

# step
step = 0.001

# evolution
evolutionTotal = []
evolutionTotal.append(totalDistance)
evolutionWorst = []
evolutionWorst.append(worstDistance)

# main loop
for i in range(number_of_iterations):

    # print iteration
    print("iteration: " + str(i+1) + "/" + str(number_of_iterations) + ", step: " + "%.6f" % round(step,6))
    
    # get a copy of the previous instrument
    newFlute = copy.deepcopy(flute)
    
    # evolve the instrument randomly
    newFlute.evolve(step)
    
    # the distance from theoretical frequency for each notes
    distances = []
    
    # for each note
    for note in fluteNotes:

        # get the playing frequency
        f1 = newFlute.playingFrequencyMinima(fingeringChartFlute[note], 1)
        f2 = newFlute.playingFrequencyMinima(fingeringChartFlute[note], 2)
        #print(note + ": " + "%.2f" % round(f1,2) + " Hz, " + "%.2f" % round(f2,2) + " Hz")
        
        newFlute.impedanceResonances[note] = (f1,f2)
        
        # theoretical frequency
        tf = theoreticalFrequency[note]
        
        if note in fluteFirstRegister:
            distances.append(abs(1200*np.log2(f1/tf)))

        elif note in fluteSecondRegister:
            distances.append(abs(1200*np.log2(f2/tf)))

        elif note in fluteThirdRegister:
            distances.append(0)
            
        # total distance
        newTotalDistance = sum(distances)

        # worst distance
        newWorstDistance = max(distances)
                        
    if newTotalDistance < totalDistance:
    
        # total distance
        totalDistance = newTotalDistance

        # worst distance
        worstDistance = newWorstDistance
        
        print("total distance from theoretical frequency: " + "%.2f" % round(totalDistance,2) + " cents")
        print("worst distance from theoretical frequency: " + "%.2f" % round(worstDistance,2) + " cents")
        flute = newFlute
        
    evolutionTotal.append(totalDistance)
    evolutionWorst.append(worstDistance)
    step -= 0.001/number_of_iterations

flute.printParameters()
flute.draw()
flute.plotReactance(fingeringChartFlute["D4"])
flute.plotImpedance(fingeringChartFlute["D4"])
flute.plotDistanceHz()
flute.plotDistanceCents()

plt.figure(figsize=(15,10))
plt.plot(evolutionTotal, "+--", color="black", label="distance entre la fréquence de l'instrument et la fréquence théorique")
plt.legend(loc="upper left")
plt.title("Évolution de la distance entre la fréquence de l'instrument et la fréquence théorique")
plt.xlabel("itérations")
plt.ylabel("distance (cents)")
plt.show()

plt.figure(figsize=(15,10))
plt.plot(evolutionWorst, "+--", color="black", label="distance de la plus mauvaise note par rapport à la fréquence théorique")
plt.legend(loc="upper left")
plt.title("Évolution de la distance de la plus mauvaise note par rapport à la fréquence théorique")
plt.xlabel("itérations")
plt.ylabel("distance (cents)")
plt.show()
