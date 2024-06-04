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

# 创建主窗口
root = tk.Tk()
root.geometry('1500x800+900+50')
root.title('势垒贯穿和隧道效应')

# 定义函数
def update():
    """
    更新界面，隐藏不需要的控件
    """
    global animation, ani
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
    """
    选择平面直角坐标系
    """
    update()
    global coordinate
    coordinate = 'cartesian'
    label1_in_frame2.config(text=f"您选择的坐标系为：平面直角坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    labelwithbutton5_in_frame1.grid(row=8, column=6)
    button5_in_frame1.grid(row=9, column=6)

def chosebutton2():
    """
    选择平面极坐标系
    """
    update()
    global coordinate
    coordinate = 'polar'
    label1_in_frame2.config(text=f"您选择的坐标系为：平面极坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    label5_in_frame1.config(text='请输入m（磁量子数）')
    label5_in_frame1.grid(row=8, column=6)
    entry4_in_frame1.grid(row=9, column=6, columnspan=3)
    label5_in_frame2.config(text='请输入m')
    label5_in_frame2.pack()

def chosebutton3():
    """
    选择柱坐标系
    """
    update()
    global coordinate
    coordinate = 'cylindrical'
    label1_in_frame2.config(text="您选择的坐标系为：柱坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    label5_in_frame1.config(text='请输入m（磁量子数）')
    label5_in_frame1.grid(row=8, column=6)
    entry4_in_frame1.grid(row=9, column=6, columnspan=3)
    label5_in_frame2.config(text='请输入m')
    label5_in_frame2.pack()
    label6_in_frame1.grid(row=10, column=6)
    entry5_in_frame1.grid(row=11, column=6, columnspan=3)
    label6_in_frame2.pack()

def chosebutton4():
    """
    选择球坐标系
    """
    update()
    global coordinate
    coordinate = 'spherical'
    label1_in_frame2.config(text="您选择的坐标系为：球坐标系。")
    label2_in_frame2.pack()
    label3_in_frame2.pack()
    label4_in_frame2.pack()
    label5_in_frame1.config(text='请输入m（磁量子数）')
    label5_in_frame1.grid(row=8, column=6)
    entry4_in_frame1.grid(row=9, column=6, columnspan=3)
    label5_in_frame2.config(text='请输入m')
    label5_in_frame2.pack()
    label7_in_frame1.grid(row=10, column=6)
    entry6_in_frame1.grid(row=11, column=6, columnspan=3)
    label7_in_frame2.pack()

def chosebutton5():
    """
    选择生成动画
    """
    global animation
    animation = True
# def chosebutton6():
#    label2_in_frame2.config(text='您选择的势垒厚度为：2a')
# def chosebutton7():
#     label2_in_frame2.config(text='您选择的势垒厚度为：3a')

def generatebutton8():
    """
    生成图像按钮的回调函数，根据选择的坐标系和参数生成图像
    """
    a = float(entry1_in_frame1.get())  # 获取势垒厚度
    label2_in_frame2.config(text='您输入的势垒厚度为：{}'.format(a))
    ratio = float(entry2_in_frame1.get())  # 获取 E/V0
    label3_in_frame2.config(text='您输入的E/V0为：{}'.format(ratio))
    V0 = float(entry3_in_frame1.get())  # 获取 V0
    label4_in_frame2.config(text='您输入的V0为：{}'.format(V0))
    m = float(entry4_in_frame1.get())  # 获取磁量子数 m

    h = float(entry5_in_frame1.get())  # 获取圆柱高度 h
    label6_in_frame2.config(text='您输入的h为：{}'.format(h))
    l = float(entry6_in_frame1.get())  # 获取角动量量子数 l
    label7_in_frame2.config(text='您输入的l为：{}'.format(l))
    global image

    if coordinate == 'cartesian':
        fig = car_plot(V0, ratio, a)  # 生成平面直角坐标系图像
        if animation:  # 如果选择了生成动画
            if not os.path.exists('image'):
                os.makedirs('image')
            if (os.path.isfile('image/gif.gif')):
                os.remove('image/gif.gif')
                print('Delete successfully')
            animation_box(V0, ratio, a)  # 生成动画
    elif coordinate == 'polar':
        label5_in_frame2.config(text='您输入的m为：{}'.format(m))
        fig = polar_plot(V0, ratio, a, m)  # 生成平面极坐标系图像
    elif coordinate == 'cylindrical':
        label5_in_frame2.config(text='您输入的m为：{}'.format(m))
        fig = cylind_plot(V0, ratio, a, h, m)  # 生成柱坐标系图像
    elif coordinate == 'spherical':
        label5_in_frame2.config(text='您输入的m为：{}'.format(m))
        if m > l:
            return messagebox2()  # 如果 m>l，弹出提示框
            fig = None
        else:
            fig = spher_plot(V0, ratio, a, m, l)  # 生成球坐标系图像
    else:
        return messagebox1()  # 如果未选择坐标系，弹出提示框

    image.get_tk_widget().destroy()  # 销毁旧图像画布
    image = FigureCanvasTkAgg(fig, master=frame3)  # 创建新图像画布
    image.draw()
    image.get_tk_widget().pack()  # 显示新图像

def animation_box(V0, ratio, a):
    """
    生成动画的函数
    """
    ani = ani_car_plot(V0, ratio, a)  # 生成动画对象
    ani.save('image/gif.gif', writer='pillow', fps=30)  # 保存动画为 gif 文件
    img = Image.open('image/gif.gif')
    def showgif(img):
        """
        显示 gif 动画的函数
        """
        top = tk.Toplevel()  # 创建顶层窗口
        top.title('一维势垒下波函数随时间变化')
        num_frames = img.n_frames  # 获取 gif 动画帧数
        frame_index = 0  # 初始化帧索引
        def animate_gif():
            """
            动画循环函数
            """
            nonlocal frame_index
            img.seek(frame_index)  # 跳转到指定帧
            photo = ImageTk.PhotoImage(img)  # 创建 PhotoImage 对象
            giflabel.config(image=photo)  # 更新标签图像
            giflabel.image = photo  # 保留对 PhotoImage 对象的引用
            frame_index = (frame_index + 1) % num_frames  # 更新帧索引
            top.after(50, animate_gif)  # 50 毫秒后再次调用动画循环函数

        giflabel = tk.Label(top)  # 创建标签用于显示动画
        giflabel.pack()
        close_button = tk.Button(top, text='关闭', command=top.destroy)  # 创建关闭按钮
        close_button.pack()
        animate_gif()  # 开始动画循环

    showgif(img)


# 定义对话框函数
def messagebox1():
    """
    未选择坐标系时的提示框
    """
    messagebox.showinfo('Note!', 'Warning! You have to choose a coordinate before generating!')

def messagebox2():
    """
    m>l 时的提示框
    """
    messagebox.showinfo('m>l')

# 创建框架
frame1 = tk.Frame(root, relief='ridge', width=400)  # frame1 用于输入参数
frame2 = tk.Frame(root, relief='ridge', width=400)  # frame2 用于显示参数信息
frame3 = tk.Frame(root, relief='ridge', width=400)  # frame3 用于显示图像
frame1.pack(side='left', fill='both')
frame2.pack(side='left', fill='both')
frame3.pack(side='left', fill='both')

# 创建控件
# frame1
label1_in_frame1 = tk.Label(frame1, text='请选择坐标系')
button1_in_frame1 = tk.Button(frame1, text='平面直角坐标', command=chosebutton1)
button2_in_frame1 = tk.Button(frame1, text='平面极坐标', command=chosebutton2)
button3_in_frame1 = tk.Button(frame1, text='柱坐标', command=chosebutton3)
button4_in_frame1 = tk.Button(frame1, text='球坐标', command=chosebutton4)
label2_in_frame1 = tk.Label(frame1, text='请输入势垒的厚度')
entry1_in_frame1 = tk.Entry(frame1, text='a')
button5_in_frame1 = tk.Checkbutton(frame1, text='动画', command=chosebutton5)
labelwithbutton5_in_frame1 = tk.Label(frame1, text='可选择生成动画')
# button6_in_frame1 = tk.Button(frame1,text='2a',command=chosebutton6)
# button7_in_frame1 = tk.Button(frame1,text='3a',command=chosebutton7)
button8_in_frame1 = tk.Button(frame1, text='生成', command=generatebutton8)
label3_in_frame1 = tk.Label(frame1, text='请输入E/V0')
entry2_in_frame1 = tk.Entry(frame1, text='E/V0')
label4_in_frame1 = tk.Label(frame1, text='请输入V0')
entry3_in_frame1 = tk.Entry(frame1, text='V0')
label5_in_frame1 = tk.Label(frame1, text='请输入m(贝塞尔函数阶数)')
entry4_in_frame1 = tk.Entry(frame1, text='m')
label6_in_frame1 = tk.Label(frame1, text='请输入h（圆柱区域的高度）')
entry5_in_frame1 = tk.Entry(frame1, text='h')
label7_in_frame1 = tk.Label(frame1, text='请输入l（角动量量子数）')
entry6_in_frame1 = tk.Entry(frame1, text='l')

# frame1.grid
label1_in_frame1.grid(row=0, column=6)
button1_in_frame1.grid(row=1, column=0, columnspan=3)
button2_in_frame1.grid(row=1, column=4, columnspan=3)
button3_in_frame1.grid(row=1, column=8, columnspan=3)
button4_in_frame1.grid(row=1, column=11, columnspan=3)
label2_in_frame1.grid(row=2, column=6)
entry1_in_frame1.grid(row=3, column=6, columnspan=3)
# button5_in_frame1.grid(row=3,column=0,columnspan=3)
# button6_in_frame1.grid(row=3,column=4,columnspan=3)
# button7_in_frame1.grid(row=3,column=8,columnspan=3)
label3_in_frame1.grid(row=4, column=6)
entry2_in_frame1.grid(row=5, column=6, columnspan=3)
label4_in_frame1.grid(row=6, column=6)
entry3_in_frame1.grid(row=7, column=6, columnspan=3)
button8_in_frame1.grid(row=15, column=6, columnspan=3)

# frame2
label1_in_frame2 = tk.Label(frame2, text='等待选择')
label2_in_frame2 = tk.Label(frame2, text='请输入势垒的厚度a')
label3_in_frame2 = tk.Label(frame2, text='请输入E/V0')
label4_in_frame2 = tk.Label(frame2, text='请输入V0')
label5_in_frame2 = tk.Label(frame2, text='请输入m')
label6_in_frame2 = tk.Label(frame2, text='请输入h')
label7_in_frame2 = tk.Label(frame2, text='请输入l')

# frame2.grid/pack
label1_in_frame2.pack()

# frame3
label1_in_frame3 = tk.Label(frame3, text='绘图展示')
label1_in_frame3.pack()

# 初始化输入框
entry1_in_frame1.insert(0, '0.2')  # a
entry2_in_frame1.insert(0, '0.2')  # ratio
entry3_in_frame1.insert(0, '1')  # V0
entry4_in_frame1.insert(0, '1')  # m
entry5_in_frame1.insert(0, '1')  # h
entry6_in_frame1.insert(0, '1')  # l

# 定义全局变量
coordinate = 0  # 初始化坐标系选择
animation = False  # 初始化动画选项
fig = plt.figure()  # 创建初始 matplotlib 图形
fig.patch.set_facecolor('white')  # 设置初始图形背景色为白色
image = FigureCanvasTkAgg(fig, master=frame3)  # 创建初始图形画布
image.draw()
image.get_tk_widget().pack()

root.mainloop()  # 进入消息循环