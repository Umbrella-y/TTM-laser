import os
import shutil
import time
#用于保存密度信息和温度信息，并加以时间戳
def save_file_with_timestamp(file_path,stamp):
    # Get the current time as a string
    #timestamp = time.strftime("%Y%m%d-%H%M%S")

    # Create a new file name with the timestamp
    file_name, file_ext = os.path.splitext(file_path)
    new_file_name = f"{file_name}-{stamp}{file_ext}"

    # Copy the file to the new name
    shutil.copy(file_path, new_file_name)

    # Return the new file path
    return new_file_name