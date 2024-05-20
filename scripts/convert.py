import os
import asyncio
from PIL import Image


QUALITY = 100
# slim later


# def get_png_files(path='.'):
#     files = []
#     for filename in os.listdir(path):
#         if filename.endswith('.png'):
#             files.append(os.path.join(path, filename))
#     return files


def recursive_get_png_files(path='.'):
    files = []
    for root, _, filenames in os.walk(path):
        for filename in filenames:
            if filename.endswith('.png'):
                files.append(os.path.join(root, filename))
    return files


async def convert(file):
    image = Image.open(file)
    image = image.convert('RGB')
    new_file = file.replace('.png', '.jpg')
    image.save(
        new_file,
        'JPEG', quality=QUALITY
    )

    print('Converted {} to {}'.format(file, new_file))

    # sleep for 10 seconds before removing the file
    await asyncio.sleep(10)
    return os.remove(file)


async def runner(tasks):
    return await asyncio.gather(*tasks)


# Convert png to jpg and /or slim it
if __name__ == '__main__':
    png_files = recursive_get_png_files('..')
    if not png_files:
        exit(0)
    convert_tasks = []
    for png_file in png_files:
        convert_tasks.append(convert(png_file))
    # loop = asyncio.get_event_loop()
    # loop.run_until_complete(asyncio.wait(tasks))
    # loop.close()
    asyncio.run(runner(convert_tasks))
