import tkinter as tk



class DrawingApp:
    def __init__(self, root):
        self.root = root
        self.canvas = tk.Canvas(root, width=400, height=400)
        self.canvas.pack()
        self.prev_x, self.prev_y = None, None
        self.coordinates = []  # To store coordinates

        # Create a "Clear" button
        self.clear_button = tk.Button(root, text="Clear", command=self.clear_canvas)
        self.clear_button.pack()

        # Create a "Get Coordinates" button
        self.get_coordinates_button = tk.Button(root, text="Get Coordinates", command=self.get_coordinates)
        self.get_coordinates_button.pack()

        # Create a "Save Coordinates" button
        self.save_coordinates_button = tk.Button(root, text="Save Coordinates", command=self.save_coordinates)
        self.save_coordinates_button.pack()

        # Draw scale labels
        self.canvas.create_text(10, 400, text="3480", anchor="sw")
        self.canvas.create_text(10, 0, text="6370", anchor="nw")

        # Bind mouse events
        self.canvas.bind("<Button-1>", self.start_drawing)
        self.canvas.bind("<B1-Motion>", self.draw)

    def start_drawing(self, event):
        self.prev_x, self.prev_y = event.x, event.y

    def draw(self, event):
        if self.prev_x is not None and self.prev_y is not None:
            x, y = event.x, event.y
            self.coordinates.append((x, y))  # Store coordinates
            self.canvas.create_line(self.prev_x, self.prev_y, x, y)
            self.prev_x, self.prev_y = x, y

    def clear_canvas(self):
        self.canvas.delete("all")  # Clear the entire canvas
        self.coordinates = []  # Clear stored coordinates

    def get_coordinates(self):
        print("Coordinates:")
        for coord in self.coordinates:
            print(f"X: {coord[0]}, Y: {coord[1]}")

    def save_coordinates(self):
        # Save coordinates to a text file with two columns (X and Y)
        with open("proto_visc/const_visc.vis", "w") as file:
            for coord in self.coordinates:
                file.write(f"{coord[0]}\n")

if __name__ == "__main__":
    root = tk.Tk()
    app = DrawingApp(root)
    root.mainloop()
