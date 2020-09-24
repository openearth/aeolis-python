import cv2
import os

fps = 5


image_folder = 'C:\\Users\\Simulations\\Blowout607'
video_name = 'C:\\Users\\Simulations\\Blowout607\\607a.avi'


images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
images.sort()
frame = cv2.imread(os.path.join(image_folder, images[0]))

height, width, layers = frame.shape

# video = cv2.VideoWriter(video_name, 0, fps, (width, height))
# video = cv2.VideoWriter(video_name, 0x00000021, fps, (width, height))
video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc('M', 'P', '4', '3'), fps, (width, height))
# video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc('M', 'J', 'P', 'G'), fps, (width, height))



for image in images:
    print(image)
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
