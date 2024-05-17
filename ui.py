from tkinter import *
from tkinter import messagebox
import os
import binascii
from Bio.Seq import Seq
from Crypto.Cipher import Salsa20
from Crypto.Random import get_random_bytes
import ipfshttpclient
import timeit
from PIL import ImageTk, Image
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from tkinter import ttk


start_widget = Tk()


start_widget.title("Data Guard")
start_widget.geometry("1300x500")

start_widget.iconbitmap('desktop_icon.ico')


start_widget.minsize(1300,500)
start_widget.maxsize(1300,500)

# Function to animate the label
def animate_label(label, x, y, direction):
    speed = 1.5  # Adjust the speed as needed
    if direction == 'right':
        x += speed
        if x >= 300:  # Adjust the center position based on your requirements
            return
    elif direction == 'left':
        x -= speed
        if x <= 200:  # Adjust the center position based on your requirements
            return
    label.place(x=x, y=y)
    start_widget.after(10, animate_label, label, x, y, direction)

# Function to start the animation
def start_animation():
    animate_label(tit_le, 200, 100, 'right')


def dna_format_conversion(_get_binary):
    apply_dna_map = {'00': 'A', '01': 'T', '10': 'G', '11': 'C'}
    
    # Pad the binary string with '0' if its length is odd
    if len(_get_binary) % 2 != 0:
        _get_binary += '0'
    
    finalized_dna_sequences = ''.join([apply_dna_map[_get_binary[i:i+2]] for i in range(0, len(_get_binary), 2)])
    return finalized_dna_sequences


def apply_dna_encryption(key):
    get_binary = bin(int.from_bytes(key.encode(), 'big'))[2:]
    finalized_dna_sequences = dna_format_conversion(get_binary)
    return str(Seq(finalized_dna_sequences).reverse_complement())


def salsa20_encrypt(key, my_plaintext_):
    key = key[:32]  # Salsa20 key size is 256 bits (32 bytes)
    cipher = Salsa20.new(key=key)
    ciphertext_ = cipher.nonce + cipher.encrypt(my_plaintext_.encode())
    return ciphertext_



def Dec(entry2,field_char):
    normal_key=field_char.get()
    ipfs_hash_=entry2.get()

    dna_encrypted_key=apply_dna_encryption(normal_key)

    print("\nDNA Encrypted Key : ",dna_encrypted_key)

    # salsa20 Encryption with DNA-encrypted key
    salsa20_key = dna_encrypted_key.encode() 

    # Retrieve the encrypted text from IPFS using the hash
    retrieved_encrypted_text, ipfs_retrieval_time = eval_(ipfs_retrieval, ipfs_hash_)

    # Perform salsa20 decryption using the original DNA key
    decrypted_text, salsa20_decryption_time = eval_(salsa20_decrypt, salsa20_key, retrieved_encrypted_text)

    print("\nDecrypted Text:", decrypted_text)
    text2.delete(1.0, END)  
    text2.insert(END, decrypted_text)

    messagebox.showinfo("Success", "Decryption successful!")

def salsa20_decrypt(key, encrypted_text):
    key = key[:32]  # Salsa20 key size is 256 bits (32 bytes)
    msg_nonce = encrypted_text[:8]
    ciphertext = encrypted_text[8:]
    cipher = Salsa20.new(key=key, nonce=msg_nonce)
    decrypted_bytes = cipher.decrypt(ciphertext)
    return decrypted_bytes

def ipfs_storage(encrypted_text):
    with ipfshttpclient.connect() as client:
        res = client.add_bytes(encrypted_text)
        return res


def eval_(func, *args, **kwargs):
    start_time = timeit.default_timer()
    result = func(*args, **kwargs)
    end_time = timeit.default_timer()
    execution_time = end_time - start_time
    return result, execution_time

def ipfs_retrieval(ipfs_hash):
    with ipfshttpclient.connect() as client:
        encrypted_text = client.cat(ipfs_hash)
        return encrypted_text


def check_existence(ipfs_hash, key):
    try:
        with open("hashes.txt", "r") as file:
            for line in file:
                saved_hash, saved_key = line.strip().split(',')
                if ipfs_hash == saved_hash and key == saved_key:
                    return True
    except FileNotFoundError:
        return False
    return False


def file_update(ipfs_hash, key):
    if not check_existence(ipfs_hash, key):
        try:
            with open("hashes.txt", "r") as file:
                get_content = file.read()
        except FileNotFoundError:
            get_content = ""

        with open("hashes.txt", "w") as file:
            file.write(f"{ipfs_hash},{key}\n")
            file.write(get_content)

def Enc(t_area1,field_char):
    try:

        normal_key=field_char.get()
        plaintext=t_area1.get("1.0",'end')

        dna_encrypted_key=apply_dna_encryption(normal_key)

        print("\nDNA Encrypted Key : ",dna_encrypted_key)

        #salsa20 Encryption with DNA-encrypted key
        salsa20_key = dna_encrypted_key.encode()  


        encrypted_text, salsa20_encryption_time = eval_(salsa20_encrypt, salsa20_key, plaintext)

        #IPFS Storage
        ipfs_hash, ipfs_storage_time = eval_(ipfs_storage, encrypted_text)


        print("\nIPFS Hash:", ipfs_hash)

        file_update(ipfs_hash, normal_key)

        messagebox.showinfo("Success", "Upload successful!")

        field_char.set("")  # Clear entry field value
        t_area1.delete("1.0", END)  # Clear text field value

    except:
        messagebox.showerror("Error", "Provide Key with correct length!")


def Go_to_home():
    global s_frame
    s_frame.pack_forget()
    s_frame = Frame(start_widget, bg="salmon")
    s_frame.pack(side="top", fill="both", expand=True)
    initial_pic = Image.open("wallpaper.jpg")
    initial_image = ImageTk.PhotoImage(initial_pic.resize((1300,500), Image.ANTIALIAS))
    start_lb = Label(s_frame, image=initial_image)
    start_lb.image = initial_image
    start_lb.pack()


    # Create a themed style for a modern look
    style = ttk.Style()
    style.configure('Title.TLabel', font=('Arial', 26, 'bold'), background='white')

    # Create the label with the themed style
    tit_le = ttk.Label(start_widget, text="Safeguarding Data: Your Trusted Data Guard.", style='Title.TLabel')
    tit_le.place(x=200, y=100)

    # Start the animation automatically when the GUI opens
    animate_label(tit_le, 200, 100, 'right')

 
def Test_Page():
    global s_frame
    s_frame.pack_forget()

    s_frame = Frame(start_widget, bg="white")
    s_frame.pack(side="top", fill="both", expand=True)

    global s1_frame
    s1_frame = Frame(s_frame, bg="#b3d9ff")
    s1_frame.place(x=0, y=0, width=650, height=500)
    s1_frame.config()

    s2_title = Label(s1_frame, text="Encryption", font="arial 16 bold", bg="#b3d9ff")
    s2_title.pack(padx=0, pady=10)

    s2_title_ = Label(s1_frame, text="Provide Key :",
                   font="arial 12 bold", bg="#b3d9ff")
    s2_title_.place(x=100, y=50)
    my_label = Label(s1_frame, text="Provide Data:", font="arial 12 bold", bg="#b3d9ff")
    my_label.place(x=100, y=120)
    
   
    global field_char,input_filed1
    field_char=StringVar()
    t_area1=Text(s1_frame,height=10,width=50)
    t_area1.place(x=100,y=150)
    input_filed1 = Entry(s1_frame, textvariable=field_char, bd=2, width=35)
    input_filed1.place(x=100, y=80)
 
    btn_up = Button(
        s1_frame, text="Encrypt",width=20,height=2, command=lambda: Enc(t_area1,field_char), bg="pink")
    btn_up.place(x=250,y=340)

    #########################################

    global widget_frame
    widget_frame = Frame(s_frame, bg="#ffe6ff")
    widget_frame.place(x=650, y=0, width=650, height=500)
    widget_frame.config()

    s2_title = Label(widget_frame, text="Decryption", font="arial 16 bold", bg="#ffe6ff")
    s2_title.pack(padx=0, pady=10)

    my_label = Label(widget_frame, text="Provide Hash:", font="arial 12 bold", bg="#ffe6ff")
    my_label.place(x=140, y=100)
    s2_title_ = Label(widget_frame, text="Provide Key :",
                   font="arial 12 bold", bg="#ffe6ff")
    s2_title_.place(x=140, y=160)

    my_label = Label(widget_frame, text="Decrypted Data:", font="arial 12 bold", bg="#ffe6ff")
    my_label.place(x=140, y=230)
   

    global v1,field_char1,text2
    field_char1=StringVar()
    v1=StringVar()

    entry2 = Entry(widget_frame, textvariable=v1, bd=2, width=25)
    entry2.place(x=280, y=100)
    input_filed1 = Entry(widget_frame, textvariable=field_char1, bd=2, width=25)
    input_filed1.place(x=280, y=160)

    text2=Text(widget_frame,height=10,width=50)
    text2.place(x=140,y=270)
 

    btn_up = Button(
        widget_frame, text="Decrypt",width=20,height=2, command=lambda: Dec(entry2,field_char1), bg="light blue")
    btn_up.place(x=460,y=140)
  

s_frame = Frame(start_widget, bg="salmon")
s_frame.pack(side="top", fill="both", expand=True)
initial_pic = Image.open("wallpaper.jpg")
initial_image = ImageTk.PhotoImage(initial_pic.resize((1300,500), Image.ANTIALIAS))
start_lb = Label(s_frame, image=initial_image)
start_lb.image = initial_image
start_lb.pack()

# Create a themed style for a modern look
style = ttk.Style()
style.configure('Title.TLabel', font=('Arial', 26, 'bold'), background='white')

# Create the label with the themed style
tit_le = ttk.Label(start_widget, text="Safeguarding Data: Your Trusted Data Guard.", style='Title.TLabel')
tit_le.place(x=200, y=100)

# Start the animation automatically when the GUI opens
animate_label(tit_le, 200, 100, 'right')


tkinter_menu = Menu(start_widget)
checkmenu = Menu(tkinter_menu)
tkinter_menu.add_command(label="My Home", command=Go_to_home)
tkinter_menu.add_command(label="Test", command=Test_Page)
start_widget.config(menu=tkinter_menu)
start_widget.mainloop()
