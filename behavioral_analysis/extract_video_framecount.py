import cv2

def get_frame_count_exact(video_path):
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        raise IOError(f"Cannot open video file: {video_path}")
    
    count = 0
    while True:
        ret, frame = cap.read()
        if not ret:
            break
        count += 1
    
    cap.release()
    return count

frame_count = get_frame_count_exact(r"E:\odorTask4Behavior\testCamera\2026_09_10_cameraTest\2026_09_10_cameraTest_g0\2026_09_10_cameraTest.mp4")
print(f"Frame count: {frame_count}")