'''
项目目标：可视化“势垒贯穿和隧道效应”
开发gui界面
展示波函数势垒贯穿和隧道效应，二维、三维图
gui界面的效果：
坐标系——直角坐标、极坐标、三维柱坐标、球坐标
四个按钮，每按一个会显示出对应的板块，可以更改势垒的边界从而改变波函数
'''
import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from plot import car_plot,polar_plot,cylind_plot,spher_plot,ani_car_plot
import os
from PIL import Image, ImageTk
#Create main window
root = tk.Tk()
root.geometry('1200x500+900+50')
root.title('势垒贯穿和隧道效应')

#Defining of Functions
#Button Widget



def update():
    global animation,ani
    label5_in_frame1.grid_forget()
    entry4_in_frame1.grid_forget()
    label5_in_frame2.pack_forget()
    label6_in_frame1.grid_forget()
    label6_in_frame2.pack_forget()
    entry5_in_frame1.grid_forget()
    entry6_in_frame1.grid_forget()
    label7_in_frame1.grid_forget()
    label7_in_frame2.pack_forget()
    labelwithbutton5_in_frame1.grid_forget()
    button5_in_frame1.grid_forget()
    animation = False
    

def chosebutton1():
    update()
    global coordinate
    coordinate = 'cartesian'
    label1_in_frame2.config(text=f"您选择的坐标系为：平面直角坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    labelwithbutton5_in_frame1.grid(row=8,column=6)
    button5_in_frame1.grid(row=9,column=6)

def chosebutton2():
    update()
    global coordinate
    coordinate = 'polar'
    label1_in_frame2.config(text=f"您选择的坐标系为：平面极坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    label5_in_frame1.config(text='请输入l（角动量量子数）')
    label5_in_frame1.grid(row=8,column=6)
    entry4_in_frame1.grid(row=9,column=6,columnspan=3)
    label5_in_frame2.config(text='请输入l')
    label5_in_frame2.pack()

def chosebutton3():
    update()
    global coordinate
    coordinate = 'cylindrical'
    label1_in_frame2.config(text="您选择的坐标系为：柱坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    label5_in_frame1.config(text='请输入l（角动量量子数）')
    label5_in_frame1.grid(row=8,column=6)
    entry4_in_frame1.grid(row=9,column=6,columnspan=3)
    label5_in_frame2.config(text='请输入l')
    label5_in_frame2.pack()
    label6_in_frame1.grid(row=10,column=6)
    entry5_in_frame1.grid(row=11,column=6,columnspan=3)
    label6_in_frame2.pack()

def chosebutton4():
    update()
    global coordinate
    coordinate = 'spherical'
    label1_in_frame2.config(text="您选择的坐标系为：球坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    label5_in_frame1.config(text='请输入m（磁量子数）')
    label5_in_frame1.grid(row=8,column=6)
    entry4_in_frame1.grid(row=9,column=6,columnspan=3)
    label5_in_frame2.config(text='请输入m')
    label5_in_frame2.pack()
    label7_in_frame1.grid(row=10,column=6)
    entry6_in_frame1.grid(row=11,column=6,columnspan=3)
    label7_in_frame2.pack()

def chosebutton5():
    global animation
    animation = True
# def chosebutton6():
#    label2_in_frame2.config(text='您选择的势垒厚度为：2a')
# def chosebutton7():
#     label2_in_frame2.config(text='您选择的势垒厚度为：3a')
def generatebutton8():
    a = entry1_in_frame1.get()
    label2_in_frame2.config(text='您输入的势垒厚度为：{}'.format(a))
    ratio = entry2_in_frame1.get()
    label3_in_frame2.config(text='您输入的E/V0为：{}'.format(ratio))
    V0 = entry3_in_frame1.get()
    label4_in_frame2.config(text='您输入的V0为：{}'.format(V0))
    m = entry4_in_frame1.get()
    
    h = entry5_in_frame1.get()
    label6_in_frame2.config(text='您输入的h为：{}'.format(h))
    l = entry6_in_frame1.get()
    label7_in_frame2.config(text='您输入的l为：{}'.format(l))
    global image
     
    if coordinate == 'cartesian':
        fig = car_plot(V0,ratio,a)
        if animation:
            if not os.path.exists('image'):
                os.makedirs('image')
            if(os.path.isfile('image/gif.gif')):
                os.remove('image/gif.gif')
                print('Delete successfully')
            animation_box(V0,ratio,a)
    elif coordinate == 'polar':
        label5_in_frame2.config(text='您输入的l为：{}'.format (m))
        fig = polar_plot(V0,ratio,a,m)
    elif coordinate == 'cylindrical':
        label5_in_frame2.config(text='您输入的l为：{}'.format (m))
        fig = cylind_plot(V0,ratio,a,h,m)
    elif coordinate == 'spherical':
        label5_in_frame2.config(text='您输入的m为：{}'.format (m))
        m,l = float(m),float(l)
        if m>l:
            return messagebox2()
            fig = None
        else:
            fig = spher_plot(V0,ratio,a,m,l)
    else:
        return messagebox1()
    
    image.get_tk_widget().destroy()
    image = FigureCanvasTkAgg(fig, master = frame3)
    image.draw()
    image.get_tk_widget().pack()

def animation_box(V0,ratio,a):
        ani = ani_car_plot(V0,ratio,a)
        ani.save('image/gif.gif',writer='pillow',fps=30)
        img = Image.open('image/gif.gif')
        def showgif(img):
            top = tk.Toplevel()
            top.title('一维势垒下波函数随时间变化')
            num_frames = img.n_frames
            frame_index = 0
            def animate_gif():
                nonlocal frame_index
                img.seek(frame_index)
                photo = ImageTk.PhotoImage(img)
                giflabel.config(image=photo)
                giflabel.image = photo
                frame_index = (frame_index+1)%num_frames
                top.after(50,animate_gif)
            giflabel = tk.Label(top)
            giflabel.pack()
            close_button = tk.Button(top,text='关闭',command=top.destroy)
            close_button.pack()
            animate_gif()
        showgif(img)


#Dialog Box
def messagebox1():
    messagebox.showinfo('Note!','Warning! You have to choose a coordinate before generating!')
def messagebox2():
    messagebox.showinfo('m>l')

#Create Frames
frame1 = tk.Frame(root,relief='ridge',width=400)#frame1用于输入
frame2 = tk.Frame(root,relief='ridge',width=400)#frame2用于显示label
frame3 = tk.Frame(root,relief='ridge',width=400)#frame3用于显示图
frame1.pack(side='left',fill='both')
frame2.pack(side='left',fill='both')
frame3.pack(side='left',fill='both')
#Create widges inside the frames
#frame1
label1_in_frame1 = tk.Label(frame1, text='请选择坐标系')
button1_in_frame1 = tk.Button(frame1,text='平面直角坐标',command=chosebutton1)
button2_in_frame1 = tk.Button(frame1,text='平面极坐标',command=chosebutton2)
button3_in_frame1 = tk.Button(frame1,text='柱坐标',command=chosebutton3)
button4_in_frame1 = tk.Button(frame1,text='球坐标',command=chosebutton4)
label2_in_frame1 = tk.Label(frame1,text='请输入势垒的厚度')
entry1_in_frame1 = tk.Entry(frame1,text='a')
button5_in_frame1 = tk.Checkbutton(frame1,text='动画',command=chosebutton5)
labelwithbutton5_in_frame1 = tk.Label(frame1,text='可选择生成动画')
# button6_in_frame1 = tk.Button(frame1,text='2a',command=chosebutton6)
# button7_in_frame1 = tk.Button(frame1,text='3a',command=chosebutton7)
button8_in_frame1 = tk.Button(frame1,text='生成',command=generatebutton8)
label3_in_frame1 = tk.Label(frame1,text='请输入E/V0')
entry2_in_frame1 = tk.Entry(frame1,text='E/V0')
label4_in_frame1 = tk.Label(frame1,text='请输入V0')
entry3_in_frame1 = tk.Entry(frame1,text='V0')
label5_in_frame1 = tk.Label(frame1,text='请输入m(贝塞尔函数阶数)')
entry4_in_frame1 = tk.Entry(frame1,text='m')
label6_in_frame1 = tk.Label(frame1,text='请输入h（圆柱区域的高度）')
entry5_in_frame1 = tk.Entry(frame1,text='h')
label7_in_frame1 = tk.Label(frame1,text='请输入l（角动量量子数）')
entry6_in_frame1 = tk.Entry(frame1,text='l')

#frame1.grid
label1_in_frame1.grid(row=0,column=6)
button1_in_frame1.grid(row=1,column=0,columnspan=3)
button2_in_frame1.grid(row=1,column=4,columnspan=3)
button3_in_frame1.grid(row=1,column=8,columnspan=3)
button4_in_frame1.grid(row=1,column=11,columnspan=3)
label2_in_frame1.grid(row=2,column=6)
entry1_in_frame1.grid(row=3,column=6,columnspan=3)
# button5_in_frame1.grid(row=3,column=0,columnspan=3)
# button6_in_frame1.grid(row=3,column=4,columnspan=3)
# button7_in_frame1.grid(row=3,column=8,columnspan=3)
label3_in_frame1.grid(row=4,column=6)
entry2_in_frame1.grid(row=5,column=6,columnspan=3)
label4_in_frame1.grid(row=6,column=6)
entry3_in_frame1.grid(row=7,column=6,columnspan=3)
button8_in_frame1.grid(row=15,column=6,columnspan=3)

#frame2
label1_in_frame2 = tk.Label(frame2, text='等待选择')
label2_in_frame2 = tk.Label(frame2,text='请输入势垒的厚度a')
label3_in_frame2 = tk.Label(frame2,text='请输入E/V0')
label4_in_frame2 = tk.Label(frame2,text='请输入V0')
label5_in_frame2 = tk.Label(frame2,text='请输入m')
label6_in_frame2 = tk.Label(frame2,text='请输入h')
label7_in_frame2 = tk.Label(frame2,text='请输入l')

#frame2.grid/pack
label1_in_frame2.pack()

#frame3
label1_in_frame3 = tk.Label(frame3, text='绘图展示')
label1_in_frame3.pack()





#initialize
coordinate = 0
entry1_in_frame1.insert(0,'0.2')#a
entry2_in_frame1.insert(0,'0.2')#ratio
entry3_in_frame1.insert(0,'1')#V0
entry4_in_frame1.insert(0,'1')#m
entry5_in_frame1.insert(0,'1')#h
entry6_in_frame1.insert(0,'1')#l

fig = plt.figure()
fig.patch.set_facecolor('white')
image = FigureCanvasTkAgg(fig, master = frame3)
image.draw()
image.get_tk_widget().pack()
animation = False


root.mainloop()