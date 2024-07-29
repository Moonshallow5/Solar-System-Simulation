from matplotlib.animation import FuncAnimation
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

class solar():
    def __init__(self,timestep,iterations):
        self.timestep=timestep
        self.iterations=iterations
        self.G=6.674e-11
        self.time=0
        
        # This intialises the list of bodies for the planets to be in
        self.body_list=[]
        # Shows the total energy of the planet after each timestep (K.e+P.e)
        self.total_energy=[]
        
        #This initialises the list of patches for the function FuncAnimation
        self.patch_list = []
        
        # This contains all the times after each timestep
        self.times = []
        
    #function that reads input from the bodies.csv file
    def read_input(self):
        self.body_list=[]
        self.time=0
        # Opens and reads all the data for the planets from the csv file
        bodies1=open('bodies.csv','r')
        bodies=bodies1.readlines()
        
        for body in bodies:
            #Split each data in the csv file with a semicolon
            
            field = body.split(";")
            
            #In the daata that was split we are ordering it for the class called Body, the values in the file are all in its S.I units
            #0=name, 1=colour, 2=mass, 3=position, 4=velocity
            self.body_list.append(Body(field[0], field[1], field[2], field[3], field[4]))
        return self
    
    #function that updates the acceleration of the body
    def calc_acceleration(self):
        
        # Resets acceleration to 0 for next planet
        acceleration=0
        for body in self.body_list:
            for other_body in self.body_list:
                #For loop checks all the planets and make sure the planet it's looking at isn't the same
                if body != other_body:
                    #sum the acceleration generated for each planet by other bodies
                    acceleration+=-self.G*(other_body.mass/(norm(body.position-other_body.position)**3))*(body.position-other_body.position)
                    
            #set the planets new acceleration
            body.new_acceleration = acceleration
            #reset acceleration to 0 for next planet
            acceleration=0
            
            
    #This functions saves the total energy of the system and updates the postion and velocity of the planets after every timestep
    def move(self):
        
        for body in range(len(self.body_list)):
            
            #Updates the postion of the planets after a timestep following the algo given
            self.body_list[body].update_position_beeman(self.timestep)
            #Calc_acceleration already loops trough all the planets so its out of the loop
        self.calc_acceleration()
        for body in range(len(self.body_list)):
            #Updates the velocity and acceleration of the planets after a timemstep following the algo given
            self.body_list[body].update_velocity_beeman(self.timestep)
            self.body_list[body].update_acceleration()
            
        #Appends the current time to the list of times
        self.times.append(self.time)
        # Appends the total energy to the list of total energies
        self.total_energy.append(self.calc_KE()+self.calc_PE())
        # Updates the time of the system
        self.time+=self.timestep
    
    #function which updates the centre of patch of the animation
    def animate(self,i):
        #write the total energy of the system to the enegies.csv file
        np.savetxt("energies.csv", self.total_energy)
        #Updates the values of the velocity and acceleration of the planets after each timestep
        self.move()
        for i in range(len(self.patch_list)):
            #set the center of each patch in animation to the position of each body
            self.patch_list[i].center = (self.body_list[i].position[0], self.body_list[i].position[1])
        return self.patch_list
        
        
    #function which does the animation of the system
    def show_animation(self):
        fig = plt.figure()
        ax = plt.axes()
        for body in range(len(self.body_list)):
                
            #Initialises the circle for each planet for the animation and the radius of the circle would be 5e9 and copy the values from the bodies.csv file
            self.patch_list.append(plt.Circle((self.body_list[body].position[0], self.body_list[body].position[1]), radius=5e9, label=self.body_list[body].name, color=self.body_list[body].colour, animated=True))
                
        for i in range(0, len(self.patch_list)):
            ax.add_patch(self.patch_list[i])
        numFrames = self.iterations
        
        #sets axis limits to +- 1.5 times the distance of the last item in the csv so the animation can fit all the planets in it
        ax.axis('scaled')
        
        lim = 1.5*self.body_list[-1].position[0]
        ax.set_xlim(-lim, +lim)
        ax.set_ylim(-lim, +lim)
        
        
        
        #This is a algo used for animation where the animate function is repeatedly called
        self.anim = FuncAnimation(fig,self.animate,numFrames,repeat=False,interval=10,blit=True)
        #set legend font size, to make it look neater
        plt.rc('legend', fontsize=6)
        #change background colour of plot to black
        ax.set_facecolor('black')
        #set title and labels for axes
        plt.title('Solar Simulation', fontweight="bold")
        plt.xlabel('X Position')
        plt.ylabel('Y Position')
        #set position of legend to outside the plot to the right
        plt.legend(bbox_to_anchor =(1.26, 0.60))
        plt.show()
        
    #function which finds the total time it takes for the planets to complete a full cycle    
    def orbital_period(self):
        for body in (self.body_list):
            #ignore orbital period of Sun
            if body.name != 'Sun':
                #move the planet up the y axis so this would be half of the orbit
                while (body.position[1] >= 0):
                    self.move()
                #move planet down the y-axis so this would then be the full orbit 
                while (body.position[1] < 0):
                    self.move()
                #for each planet print the name and time taken for orbit
                print(str(body.name) + " orbital period = " + str(self.time) + " seconds" +" which is " +str(round((self.time/31536000),2))+" Earth years")
        
        
    #function that calculates the kinetic energy of the system
    def calc_KE(self):
        kinetic_energy = []
        for body in range (len(self.body_list)):
            #convert velocity vector to speed due to 1/2mv^2 and I used np.linalg.norm for it
            kinetic_energy.append((1/2)*self.body_list[body].mass*np.linalg.norm(self.body_list[body].velocity)**2)
        #sum the kinetic energies of each planet
        total_kinetic_energy = sum(kinetic_energy)
        return total_kinetic_energy
    
    #function that calculates the potential energy of the system
    def calc_PE(self):
        potential_energies=[]
        potential_energy=0
        for body in self.body_list:
            for other_body in self.body_list:
                if(body!=other_body):
                #convert position vector to distance due to -GmM/r and I used np.linalg.norm for it again
                    potential_energy += (self.G*body.mass*other_body.mass)/np.linalg.norm((other_body.position-body.position))
            potential_energies.append(potential_energy)
            potential_energy=0
            
        total_potential_energy = -(1/2)*sum(potential_energies)
        return total_potential_energy
    
    #function that plots the total energy of the system with the timestep
    def plot_total_energy(self):
        #line plot of time against total energy of system
        plt.plot(self.times, self.total_energy)
        #this sets title and labels for axes
        plt.xlabel("Time (seconds)")
        plt.title("Total Energy of System over Time")
        plt.ylabel("Energy (Joules)")
        plt.show()
        
        #This method is just for Experiment 3
    def satelite_to_Mars(self):
        self.times = []
        #list containing the distance of satellite mars at each timestep
        distances_satellite_mars = []
        for body in self.body_list:
            for other_body in self.body_list:
                if body.name=='Perseverance' and other_body.name=='Mars':
                    # The value below is the total time for the full orbit of Mars
                    while self.time < 59356800:
                        #convert position vector to distance using linalg.norm again
                        distances_satellite_mars.append(np.linalg.norm(body.position-other_body.position))
                        self.move()
        #print the minimum distance satellite gets to mars in kilometres
        print(str("Closest distance of satellite and Mars = ")+str((min(distances_satellite_mars)/1000).round(2))+" kilometres")
        print("Time to reach closest distance to Mars = " +str(self.times[distances_satellite_mars.index(min(distances_satellite_mars))])+ " seconds")




            
        
    
   
        
    
            
        
        
    
class Body(object):
    def __init__(self,name,colour,mass,position,velocity,new_acceleration=[0,0],current_acceleration=[0,0],previous_acceleration=[0,0]):
        self.name=name
        self.colour=colour
        self.mass=float(mass)
        # converts the data in the [] to be an np.array
        self.position=np.array(position.strip(' [  ]').split(','),dtype=float)
        
        self.velocity = np.array(velocity.strip('[ ]\n').split(','), dtype=float)
        self.new_acceleration=np.array(new_acceleration,dtype=float)
        self.current_acceleration=np.array(current_acceleration,dtype=float)
        self.previous_acceleration=np.array(previous_acceleration,dtype=float)
        
        
        
    def update_position_beeman(self,timestep):
        self.position=self.position+self.velocity*timestep+(1/6)*(4*self.current_acceleration-self.previous_acceleration)*timestep**2
    
    def update_velocity_beeman(self,timestep):
        self.velocity=self.velocity+(1/6)*((2*self.new_acceleration)+(5*self.current_acceleration)-(self.previous_acceleration))*timestep
        
        #Algo for updating position using forward euler for every timestep
    def update_position_euler(self,timestep):
        self.position=self.position+self.velocity*timestep
        
        #Algo for updating velocty using forward euler for every timestep
    def update_velocity_euler(self,timestep):
        self.velocity=self.velocity+self.current_acceleration*timestep
        
    def update_acceleration(self):
        self.previous_acceleration=self.current_acceleration
        self.current_acceleration=self.new_acceleration
    
    
    
    
def main():
    a=solar(86400,700)

    a.read_input().show_animation()
    a.read_input().plot_total_energy()
    a.read_input().orbital_period()
    
    a.read_input().satelite_to_Mars()
    
main()