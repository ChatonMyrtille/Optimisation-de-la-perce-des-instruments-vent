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
            
                # when the hole is part of the bell
                if currentToneHole == 19 or currentToneHole == 20:
                
                    # update the transfer matrix
                    T = np.dot(T, cell.transferMatrix(f, 1))
                    
                # otherwise
                else:
                
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
                
        # radiation impedance
        Zr = complex(0,0)
        T = np.dot(T,np.matrix([[Zr],[complex(1,0)]], dtype=complex))
        
        # input impedance
        Zin = T.item(0)/T.item(1)
        return Zin
        
    def draw(self):
        image = Image.new("RGBA", (1000, 500), (255, 255, 255, 255))
        draw = ImageDraw.Draw(image)
        fontSmall = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 10)
        fontLarge = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 15)
        white = (255,255,255,255)
        black = (0,0,0,255)
        
        # draw cylinders and cones
        currentLength = 100
        for i in range(len(self.cells)):
        
            # draw cylinders
            if self.cells[i].__class__.__name__ == "CylCell":
                        
                # first cylinder
                if i == 0:
                
                    # draw exterior cylinder
                    draw.line((currentLength, 250 - (self.cells[i].a+0.0075)*1000, currentLength + (self.cells[i].L)*1000, 250 - (self.cells[i].a+0.0075)*1000), fill=black)
                    
                    # draw interior cylinder
                    draw.line((currentLength, 250 - self.cells[i].a*1000, currentLength + (self.cells[i].L)*1000, 250 - self.cells[i].a*1000), fill=black)
                    
                # last cylinder
                elif i == len(self.cells)-1:
                    # draw exterior cylinder
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000, currentLength + (self.cells[i].L)*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000), fill=black)
                    
                    # draw interior cylinder
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - self.cells[i].a*1000, currentLength + (self.cells[i].L)*1000, 250 - self.cells[i].a*1000), fill=black)
                
                # cylinder - cylinder - cylinder
                elif self.cells[i-1].__class__.__name__ == "CylCell" and self.cells[i+1].__class__.__name__ == "CylCell":
                    
                    # draw exterior cylinder
                    draw.line((currentLength, 250 - (self.cells[i-1].a+0.0075)*1000, currentLength + (self.cells[i].L)*1000, 250 - (self.cells[i+1].a+0.0075)*1000), fill=black)
                    
                    # draw interior cylinder
                    draw.line((currentLength, 250 - self.cells[i-1].a*1000, currentLength + (self.cells[i].L)*1000, 250 - self.cells[i+1].a*1000), fill=black)
                
                # cylinder - cylinder - hole
                elif self.cells[i-1].__class__.__name__ == "CylCell" and self.cells[i+1].__class__.__name__ == "HoleCell":
                
                    # draw exterior left cylinder
                    draw.line((currentLength, 250 - (self.cells[i].a+0.0075)*1000, currentLength + (self.cells[i].L/2)*1000, 250 - (self.cells[i].a+0.0075)*1000), fill=black)
                    
                    # draw exterior right cylinder
                    draw.line((currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+self.cells[i+1].h)*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - (self.cells[i].a+self.cells[i+1].h)*1000), fill=black)
                    
                    # draw exterior step cylinder
                    draw.line((currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+0.0075)*1000 ,currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+self.cells[i+1].h)*1000), fill=black)
                    
                    # draw interior cylinder
                    draw.line((currentLength, 250 - self.cells[i].a*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - self.cells[i].a*1000), fill=black)
                    
                    
                # hole - cylinder - cylinder
                elif self.cells[i-1].__class__.__name__ == "HoleCell" and self.cells[i+1].__class__.__name__ == "CylCell":
                
                    # draw exterior left cylinder
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000, currentLength + (self.cells[i].L/2)*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000), fill=black)
                    
                    # draw exterior right cylinder
                    draw.line((currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+0.0075)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i].a+0.0075)*1000), fill=black)
                    
                    # draw exterior step cylinder
                    draw.line((currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000 ,currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+0.0075)*1000), fill=black)
                    
                    # draw interior cylinder
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - self.cells[i].a*1000, currentLength + self.cells[i].L*1000, 250 - self.cells[i].a*1000), fill=black)
                    
                # other cylinder
                else:
                    # draw exterior left cylinder
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000, currentLength + (self.cells[i].L/2)*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000), fill=black)
                    
                    # draw exterior right cylinder
                    draw.line((currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+self.cells[i+1].h)*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - (self.cells[i].a+self.cells[i+1].h)*1000), fill=black)
                    
                    # draw exterior step cylinder
                    draw.line((currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+self.cells[i-1].h)*1000 ,currentLength + self.cells[i].L/2*1000, 250 - (self.cells[i].a+self.cells[i+1].h)*1000), fill=black)
                    
                    # draw interior cylinder
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - self.cells[i].a*1000, currentLength + (self.cells[i].L - self.cells[i+1].b)*1000, 250 - self.cells[i].a*1000), fill=black)
                
                
                draw.line((currentLength, 250 + self.cells[i].a*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a*1000), fill=black)
                draw.line((currentLength, 250 + (self.cells[i].a+0.0075)*1000, currentLength + self.cells[i].L*1000, 250 + (self.cells[i].a+0.0075)*1000), fill=black)
                
                # increment current length
                currentLength += self.cells[i].L*1000
            
            # draw cones
            elif self.cells[i].__class__.__name__ == "ConeCell":
            
                # hole - [cone] - cone
                if self.cells[i-1].__class__.__name__ == "HoleCell" and self.cells[i+1].__class__.__name__ == "ConeCell":
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - (self.cells[i-1].a + self.cells[i-1].h)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i-1].a + self.cells[i-1].h)*1000), fill=black)
                    draw.line((currentLength + self.cells[i-1].b*1000, 250 - self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 - self.cells[i].a2*1000), fill=black)
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                    draw.line((currentLength, 250 + (2*0.0075)*1000, currentLength + self.cells[i].L*1000, 250 + (2*0.0075)*1000), fill=black)
                    currentLength += self.cells[i].L*1000
            
                # cone - [cone] - cone
                if self.cells[i-1].__class__.__name__ == "ConeCell" and self.cells[i+1].__class__.__name__ == "ConeCell":
                    # hole - cone - [cone]
                    if self.cells[i-2].__class__.__name__ == "HoleCell":
                        draw.line((currentLength, 250 - (self.cells[i-2].a + self.cells[i-2].h)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i-2].a + self.cells[i-2].h)*1000), fill=black)
                        # draw step
                        draw.line((currentLength + self.cells[i].L*1000, 250 - (self.cells[i-2].a + self.cells[i-2].h)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i+3].a + self.cells[i+3].h)*1000), fill=black)
                    # [cone] - cone - hole
                    if self.cells[i+2].__class__.__name__ == "HoleCell":
                        draw.line((currentLength, 250 - (self.cells[i+2].a + self.cells[i+2].h)*1000, currentLength + self.cells[i].L*1000, 250 - (self.cells[i+2].a + self.cells[i+2].h)*1000), fill=black)
                    draw.line((currentLength, 250 - self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 - self.cells[i].a2*1000), fill=black)
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                    draw.line((currentLength, 250 + (2*0.0075)*1000, currentLength + self.cells[i].L*1000, 250 + (2*0.0075)*1000), fill=black)
                    currentLength += self.cells[i].L*1000
                
                # cone - [cone] - hole
                if self.cells[i-1].__class__.__name__ == "ConeCell" and self.cells[i+1].__class__.__name__ == "HoleCell":
                    draw.line((currentLength, 250 - (self.cells[i+1].a + self.cells[i+1].h)*1000, currentLength + self.cells[i].L*1000 - self.cells[i+1].b*1000, 250 - (self.cells[i+1].a + self.cells[i+1].h)*1000), fill=black)
                    draw.line((currentLength, 250 - self.cells[i].a1*1000, currentLength + self.cells[i].L*1000 - self.cells[i+1].b*1000, 250 - self.cells[i].a2*1000), fill=black)
                    draw.line((currentLength, 250 + self.cells[i].a1*1000, currentLength + self.cells[i].L*1000, 250 + self.cells[i].a2*1000), fill=black)
                    draw.line((currentLength, 250 + (2*0.0075)*1000, currentLength + self.cells[i].L*1000, 250 + (2*0.0075)*1000), fill=black)
                    currentLength += self.cells[i].L*1000
        
        # draw holes
        currentLength = 100
        
        for cell in self.cells:
            if cell.__class__.__name__ == "CylCell" or cell.__class__.__name__ == "ConeCell":
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
        
        for f in range(100, 3000):
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
    def evolve(self):
        for i in range(5, len(self.cells)-5):

            if self.cells[i].__class__.__name__ == "CylCell":
            
                # change cylinder length
                self.cells[i].L += rnd.random()*0.0010 - 0.0005
                if self.cells[i].L < 0.010: self.cells[i].L = 0.010
                if self.cells[i].L > 0.050: self.cells[i].L = 0.050
                
            if self.cells[i].__class__.__name__ == "HoleCell":
            
                # change hole radius
                self.cells[i].b += rnd.random()*0.0010 - 0.0005
                if self.cells[i].b < 0.0010: self.cells[i].b = 0.0010
                if self.cells[i].b > 0.0050: self.cells[i].b = 0.0050
                
                # change chimney length
                self.cells[i].h += rnd.random()*0.0010 - 0.0005
                if self.cells[i].h < 0.0010: self.cells[i].h = 0.0010
                if self.cells[i].h > 0.0050: self.cells[i].h = 0.0050
    
    # plot distance in Hz
    def plotDistanceHz(self):
        x = ["Ré2", "Ré2#", "Mi2", "Fa2", "Fa2#", "Sol2", "Sol2#", "La2", "La2#", "Si2", "Do3", "Do3#", "Ré3", "Ré3#", "Mi3", "Fa3", "Fa3#", "Sol3", "Sol3#", "La3", "La3#", "Si3", "Do4", "Do4#", "Ré4", "Ré4#", "Mi4", "Fa4"]
        y1 = []
        y2 = []
        y3 = []
        y4 = []

        for note in clarinetNotes:
            tf = theoreticalFrequency[note]
            f1 = self.impedanceResonances[note][0]
            f2 = self.impedanceResonances[note][1]
        
            # x.append(note)
            y1.append(tf)
            y2.append(tf*3)
            y3.append(f1)
            y4.append(f2)
        
        plt.figure(figsize=(15,10))
        plt.plot(x, y1, "-", color="grey", label="fréquence théorique")
        plt.plot(x, y2, "-", color="grey")
        plt.plot(x, y3, "+--", color="black", label="maxima de l'impedance")
        plt.plot(x, y4, "+--", color="black")
        plt.legend(loc="upper left")
        plt.title("Distance entre la fréquence de l'instrument et la fréquence théorique")
        plt.xlabel("notes")
        plt.ylabel("fréquence (Hz)")
        plt.show()
    
    # plot distance in cents
    def plotDistanceCents(self):
        x = ["Ré2", "Ré2#", "Mi2", "Fa2", "Fa2#", "Sol2", "Sol2#", "La2", "La2#", "Si2", "Do3", "Do3#", "Ré3", "Ré3#", "Mi3", "Fa3", "Fa3#", "Sol3", "Sol3#", "La3", "La3#", "Si3", "Do4", "Do4#", "Ré4", "Ré4#", "Mi4", "Fa4"]
        y = []
        z = []

        for note in clarinetNotes:
            # x.append(note)
            z.append(0)
            
            f1 = self.impedanceResonances[note][0]
            f2 = self.impedanceResonances[note][1]
            
            tf = theoreticalFrequency[note]
            
            if note in clarinetFirstRegister:
                a1 = tf
                b1 = f1
                a2 = tf*3
                b2 = f2
                y.append((1200*np.log2(b1/a1) + 1200*np.log2(b2/a2)) / 2)
            elif note in clarinetSecondRegister:
                a = tf
                b = f2
                y.append(1200*np.log2(b/a))
        
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
        currentLength = 0
        holeNumber = 1
        print("hole # \t\tposition \thole diameter \tchimney length")
        for i in range(len(self.cells)):
            if self.cells[i].__class__.__name__ == "CylCell" or self.cells[i].__class__.__name__ == "ConeCell":
                currentLength += self.cells[i].L
            if self.cells[i].__class__.__name__ == "HoleCell":
                if holeNumber == 1:
                    print(str(holeNumber) + "st hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                if holeNumber == 2:
                    print(str(holeNumber) + "rd hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                if holeNumber == 19:
                    print("1st vent-hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                if holeNumber == 20:
                    print("2nd vent-hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                if holeNumber >= 3 and holeNumber <= 18:
                    print(str(holeNumber) + "th hole:\t" + "%.2f" % round(currentLength*1000,2) + " mm" + "\tø " + "%.2f" % round(self.cells[i].b*2*1000,2) + " mm" + "\t" + "%.2f" % round(self.cells[i].h*1000,2) + " mm")
                holeNumber += 1
    
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
    return np.matrix([[(a2/a1)*cmath.cos(kc*L)-cmath.sin(kc*L)/(k*x1),complex(0,Zc*cmath.sin(kc*L))],[(1.0/Zc)*(complex(0, (1+(1.0/((k**2)*x1*x2)))*cmath.sin(kc*L))+complex(0,-(1.0/x1-1.0/x2)*cmath.cos(kc*L)/k)),(a1/a2)*cmath.cos(kc*L)+cmath.sin(kc*L)+cmath.sin(kc*L)/(k*x2)]], dtype=complex)

# bore radius (7.5 mm)
a = 0.0075

# hole radius (4.5 mm)
b = 0.0045

# distance between 2 holes (20 mm)
L = 0.020

# chimney length (3.5 mm)
h = 0.0035

# create a empty instrument
clarinet = Instrument()

# bore (150 mm)
clarinet.addCell(CylCell(a, 0.100))
clarinet.addCell(CylCell(a, 0.025))
clarinet.addCell(CylCell(a, 0.025))

# register hole
clarinet.addCell(HoleCell(a, 0.0010, 0.015))

# bore (50 mm)
clarinet.addCell(CylCell(a, 0.025))
clarinet.addCell(CylCell(a, 0.025))

# 2nd hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 3rd hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 4th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 5th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 6th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 7th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 8th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 9th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 10th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 11th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 12th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 13th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 14th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 15th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 16th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 17th hole
clarinet.addCell(HoleCell(a, b, h))

# bore (20 mm)
clarinet.addCell(CylCell(a, L))

# 18th hole
clarinet.addCell(HoleCell(a, b, h))

# bell
clarinet.addCell(CylCell(a, 2*0.0182))
clarinet.addCell(HoleCell(a, 0.004, 0.010))
clarinet.addCell(CylCell(a, 2*0.0182))
clarinet.addCell(HoleCell(a, 0.004, 0.010))
clarinet.addCell(CylCell(a, 0.0182))

# fingering chart for an 18 tone holes clarinet
fingeringChartClarinet = {}
fingeringChartClarinet["D3"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
fingeringChartClarinet["D#3"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
fingeringChartClarinet["E3"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1]
fingeringChartClarinet["F3"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1]
fingeringChartClarinet["F#3"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
fingeringChartClarinet["G3"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1]
fingeringChartClarinet["G#3"] = [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1]
fingeringChartClarinet["A3"] = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1]
fingeringChartClarinet["A#3"] = [0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]
fingeringChartClarinet["B3"] = [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["C4"] = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["C#4"] = [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["D4"] = [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["D#4"] = [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["E4"] = [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["F4"] = [0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["F#4"] = [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["G4"] = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["G#4"] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
fingeringChartClarinet["A4"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
fingeringChartClarinet["A#4"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
fingeringChartClarinet["B4"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1]
fingeringChartClarinet["C5"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1]
fingeringChartClarinet["C#5"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
fingeringChartClarinet["D5"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1]
fingeringChartClarinet["D#5"] = [1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1]
fingeringChartClarinet["E5"] = [1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1]
fingeringChartClarinet["F5"] = [1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]

# theoretical frequencies
theoreticalFrequency = {}
theoreticalFrequency["D3"] = 146.8324
theoreticalFrequency["D#3"] = 155.5635
theoreticalFrequency["E3"] = 164.8138
theoreticalFrequency["F3"] = 174.6141
theoreticalFrequency["F#3"] = 184.9972
theoreticalFrequency["G3"] = 195.9977
theoreticalFrequency["G#3"] = 207.6523
theoreticalFrequency["A3"] = 220.0000
theoreticalFrequency["A#3"] = 233.0819
theoreticalFrequency["B3"] = 246.9417
theoreticalFrequency["C4"] = 261.6256
theoreticalFrequency["C#4"] = 277.1826
theoreticalFrequency["D4"] = 293.6648
theoreticalFrequency["D#4"] = 311.1270
theoreticalFrequency["E4"] = 329.6276
theoreticalFrequency["F4"] = 349.2282
theoreticalFrequency["F#4"] = 369.9944
theoreticalFrequency["G4"] = 391.9954
theoreticalFrequency["G#4"] = 415.3047
theoreticalFrequency["A4"] = 440.0000
theoreticalFrequency["A#4"] = 466.1638
theoreticalFrequency["B4"] = 493.8833
theoreticalFrequency["C5"] = 523.2511
theoreticalFrequency["C#5"] = 554.3653
theoreticalFrequency["D5"] = 587.3295
theoreticalFrequency["D#5"] = 622.2540
theoreticalFrequency["E5"] = 659.2551
theoreticalFrequency["F5"] = 698.4565
theoreticalFrequency["F#5"] = 739.9888
theoreticalFrequency["G5"] = 783.9909
theoreticalFrequency["G#5"] = 830.6094
theoreticalFrequency["A5"] = 880.0000
theoreticalFrequency["A#5"] = 932.3275
theoreticalFrequency["B5"] = 987.7666
theoreticalFrequency["C6"] = 1046.502
theoreticalFrequency["C#6"] = 1108.731
theoreticalFrequency["D6"] = 1174.659
theoreticalFrequency["D#6"] = 1244.508
theoreticalFrequency["E6"] = 1318.510
theoreticalFrequency["F6"] = 1396.913
theoreticalFrequency["F#6"] = 1479.978
theoreticalFrequency["G6"] = 1567.982
theoreticalFrequency["G#6"] = 1661.219
theoreticalFrequency["A6"] = 1760.000

# range of the clarinet
clarinetNotes = ["D3", "D#3", "E3", "F3", "F#3", "G3", "G#3", "A3", "A#3", "B3", "C4", "C#4", "D4", "D#4", "E4", "F4", "F#4", "G4", "G#4", "A4", "A#4", "B4", "C5", "C#5", "D5", "D#5", "E5", "F5"]

clarinetFirstRegister = ["D3", "D#3", "E3", "F3", "F#3", "G3", "G#3", "A3", "A#3", "B3", "C4", "C#4", "D4", "D#4", "E4", "F4", "F#4", "G4", "G#4"]

clarinetSecondRegister = ["A4", "A#4", "B4", "C5", "C#5", "D5", "D#5", "E5", "F5"]

# the total distance from theoretical frequency
totalDistance = 0

# for each note
for note in clarinetNotes:

    # playing frequency
    f1 = clarinet.playingFrequencyMaxima(fingeringChartClarinet[note], 1)
    f2 = clarinet.playingFrequencyMaxima(fingeringChartClarinet[note], 2)
    #print("playing frequency: " + "%.2f" % round(pf,2) + "Hz")
    
    clarinet.impedanceResonances[note] = (f1,f2)
    
    # theoretical frequency
    tf = theoreticalFrequency[note]
        
    # distances from theoretical frequencies
    d11 = math.sqrt((f1-tf)**2)
    d23 = math.sqrt((f2-tf*3)**2)
    d21 = math.sqrt((f2-tf)**2)
    
    if note in clarinetFirstRegister:
        totalDistance += (d11 + d23) / 2
            
    elif note in clarinetSecondRegister:
        totalDistance += d21
    
bestDistance = totalDistance
print("distance from theoretical frequency: " + "%.2f" % round(bestDistance,2) + " Hz")

# number of iterations
steps = 300

# array of distances
distances = []
distances.append(bestDistance)

# main loop
for i in range(steps):

    # print steps
    print(str(i+1) + "/" + str(steps))
    
    # get a copy of the previous instrument
    newClarinet = copy.deepcopy(clarinet)
    
    # evolve the instrument randomly
    newClarinet.evolve()
    
    # update total distance from theoretical frequency
    totalDistance = 0
    
    # for each note
    for note in clarinetNotes:

        # playing frequency
        f1 = newClarinet.playingFrequencyMaxima(fingeringChartClarinet[note], 1)
        f2 = newClarinet.playingFrequencyMaxima(fingeringChartClarinet[note], 2)
        #print("playing frequency: " + "%.2f" % round(pf,2) + "Hz")
        
        newClarinet.impedanceResonances[note] = (f1,f2)
    
        # theoretical frequency
        tf = theoreticalFrequency[note]
        
        # distances from theoretical frequencies
        d11 = math.sqrt((f1-tf)**2)
        d23 = math.sqrt((f2-tf*3)**2)
        d21 = math.sqrt((f2-tf)**2)
    
        if note in clarinetFirstRegister:
            totalDistance += (d11 + d23) / 2
            
        elif note in clarinetSecondRegister:
            totalDistance += d21
                        
    if (totalDistance < bestDistance):
        bestDistance = totalDistance
        print("distance from theoretical frequency: " + "%.2f" % round(bestDistance,2) + " Hz")
        clarinet = newClarinet
    
    distances.append(bestDistance)

clarinet.printParameters()
clarinet.draw()
clarinet.plotReactance(fingeringChartClarinet["D3"])
clarinet.plotImpedance(fingeringChartClarinet["D3"])
clarinet.plotDistanceHz()
clarinet.plotDistanceCents()


plt.figure(figsize=(15,10))
plt.plot(distances, "+--", color="black", label="distance entre la fréquence de l'instrument et la fréquence théorique")
plt.legend(loc="upper left")
plt.title("Évolution de la distance entre la fréquence de l'instrument et la fréquence théorique")
plt.xlabel("itérations")
plt.ylabel("distance (Hz)")
plt.show()

