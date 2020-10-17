import argparse
import numpy as np
import os
import shutil
import sqlite3
import types
import sys
IS_PYTHON3 = sys.version_info[0] >= 3

def array_to_blob(array):
    if IS_PYTHON3:
        return array.tostring()
    else:
        return np.getbuffer(array)

def recover_database_images_and_ids(database_path):
    # Connect to the database.
    connection = sqlite3.connect(database_path)
    cursor = connection.cursor()

    # Recover database images and ids.
    images = {}
    cameras = {}
    cursor.execute("SELECT name, image_id, camera_id FROM images;")
    for row in cursor:
        images[row[0]] = row[1]
        cameras[row[0]] = row[2]

    # Close the connection to the database.
    cursor.close()
    connection.close()

    return images, cameras


def import_line_segments(images, database_path, line_segments_path):
    # Connect to the database.
    connection = sqlite3.connect(database_path)
    cursor = connection.cursor()
    
    #cursor.execute("DELETE FROM keypoints;")
    #cursor.execute("DELETE FROM descriptors;")
    #cursor.execute("DELETE FROM matches;")
    #connection.commit()

    # Import the features.    
    f = open(line_segments_path)
    
    while True:

        line = f.readline()

        if line == '':
            break

        tmp = line.split(' ')
        assert(len(tmp) == 2)

        image_name = tmp[0]
        num_lines = int(tmp[1])

        image_id = images[image_name]

        line_segs = np.empty((num_lines,4))

        for i in range(0,num_lines):
            line = f.readline()
            tmp = line.split(' ')
            line_segs[i,0] = float(tmp[0])
            line_segs[i,1] = float(tmp[1])
            line_segs[i,2] = float(tmp[2])
            line_segs[i,3] = float(tmp[3])

        
        line_segs = line_segs.astype(np.float32)
        line_segs_str = line_segs.tostring()
        cursor.execute("DELETE FROM line_segments WHERE image_id = ?;", (image_id,))
        cursor.execute("INSERT INTO line_segments(image_id, rows, cols, data) VALUES(?, ?, ?, ?);",
                       (image_id, line_segs.shape[0], line_segs.shape[1], line_segs_str))
        connection.commit()
        
    
    # Close the connection to the database.
    cursor.close()
    connection.close()


def import_line_matches(images, database_path, matches_path):
    # Connect to the database.
    connection = sqlite3.connect(database_path)
    cursor = connection.cursor()
        
    # Import the features.    
    f = open(matches_path)
    
    while True:

        line = f.readline()

        if line == '':
            break

        tmp = line.split(' ')
        assert(len(tmp) == 3)

        image1_name = tmp[0]
        image2_name = tmp[1]
        num_matches = int(tmp[2])

        image1_id = images[image1_name]
        image2_id = images[image2_name]

        matches = np.empty((num_matches,2))

        for i in range(0,num_matches):
            line = f.readline()
            tmp = line.split(' ')
            matches[i,0] = int(tmp[0])
            matches[i,1] = int(tmp[1])
            

        if image1_id > image2_id:
            matches = matches[:, [1, 0]]
        
        image_pair_id = image_ids_to_pair_id(image1_id, image2_id)
        
        matches = matches.astype(np.uint32)        
        matches_str = matches.tostring()

        cursor.execute("DELETE FROM line_matches WHERE pair_id = ?;", (image_pair_id,))
        
        cursor.execute("INSERT INTO line_matches(pair_id, rows, cols, data) VALUES(?, ?, ?, ?);",
                       (image_pair_id, matches.shape[0], matches.shape[1], matches_str))
        connection.commit()

        
    
    # Close the connection to the database.
    cursor.close()
    connection.close()


def image_ids_to_pair_id(image_id1, image_id2):
    if image_id1 > image_id2:
        return 2147483647 * image_id2 + image_id1
    else:
        return 2147483647 * image_id1 + image_id2



if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('--database_path', required=True, help='Name of the COLMAP database')
    parser.add_argument('--line_segments', required=False, help='File containing line segments')
    parser.add_argument('--matches', required=False, help='File containing matches')
    
    parser

    args = parser.parse_args()
        
    # Reconstruction pipeline.
    images, cameras = recover_database_images_and_ids(args.database_path)

    if args.line_segments is not None:
        print('==> Importing line segments')
        import_line_segments(images, args.database_path, args.line_segments)

    if args.matches is not None:
        print('==> Importing matches...')
        import_line_matches(images, args.database_path, args.matches)
    
    print('Done!')