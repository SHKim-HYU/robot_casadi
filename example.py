from robot import Robot

robot_name1 = "ASM_CAESAR_noload"
robot_name2 = "ASM_CAESAR"

rob1 = Robot(name=robot_name1, analytical_derivatives=False)
rob2 = Robot(name=robot_name2, analytical_derivatives=False)

print(rob1.M([0,0,0,0,0,0,0]))
print(rob1.C([0,0,0,0,0,0,0],[0,0,10,0,0,0,0]))
print(rob1.G([0,1.5709,0,0,0,0,0]))

print(rob2.M([0,0,0,0,0,0,0]))
print(rob2.C([0,0,0,0,0,0,0],[0,0,10,0,0,0,0]))
print(rob2.G([0,1.5709,0,0,0,0,0]))

