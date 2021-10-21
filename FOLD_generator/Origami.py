#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 09:44:18 2020

@author: evbernardes
"""
import numpy as np
import json
import collections
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Origami:
    """ Class defining origami vertices, faces and folds in order to create 
    a FOLD file. More info on all of the possible variables of a FOLD file
    (including the ones not used in this class) at:
    https://github.com/edemaine/fold/blob/master/doc/spec.md

    Attributes
    ---------
    Metadata:
        file_version: int
                define version of FOLD file
        file_creator: string
                writes the software used to generate the file
        file_author: string
                Creator of the file
                
    FOLD config options:
        file_classes: list of strings
                = 'singleModel'
        frame_classes: list of strings
                = 'foldedForm'
        frame_attributes: list of strings
                = ['3D','orientable']
        frame_unit: string
                = 'unit'
        frame_title: string
                title of frame
    
    Vertices variables:
        vertices_coords: list of triples
            list of coordinates of each vertex
        
        vertices_names: list of strings
            (Helper variable) list of custom names for each vertex
        
        vertices_size: int
            number of vertices
            
        idx: dict
            for any given vertex name, gives the index on vertices_coords of
            the vertex which has this name. Ex:
                idx = {'square_upperLeft': 0,
                       'square_upperRight': 1,
                       'square_bottomRight': 2,
                       'square_bottomLeft': 3}
    Edges variables:
        edges_vertices: list of duplets
            each element contains the indexes of two vertices connected by an
            edge
            
        edges_length: list
            length of each edge
            
        edges_assignment: list of strings
            ith element defines the type of the ith edge in edges_vertices.
            Available options:
                'M' = mountain fold
                'V' = valley fold
                'B' = boundary
                'U' = unknown
                'F' = flat
                
        edges_foldAngle: list
            defines fixed angle of a fold
                
    Faces variables:
        faces_vertices: list of lists of ints
            each element is a list containing the indexes of vertices defining
            a single face
            
        faceOrders: list of lists  
            defines ordering of faces when necessary (which face is 
            below, etc.). For more info, see:
            https://github.com/edemaine/fold/blob/master/doc/spec.md
        
        faces_size: int
            number of faces
    

    

    Methods
    ---------
    add_vertex(name,coords)
    - name: string
    - coords: duplet
        Adds "coords" to self.vertices_coords, "name" to self.vertices_names 
        and create a new entry in self.idx using "name" to point to this new
        vertex. Also updates self.vertices_size
        
    connect(name_1,name_2,assignment,angle = 'null')
    - name_1: string
    - name_2: string
    - assignement: string
        Gets indices of vertices named "name_1" and "name_2" and create an
        edge between these two vertices in self.edges_vertices. Then,
        adds "assignment" to self.edges_assignments

    add_face(vertices_names)
    - vertices_names: list of strings
        Gets indices of vertices with names in "vertices_names" and create
        a face with these vertices, adding it to self.faces_vertices. Updates
        self.faces_size
        
    plot(lines=True, faces = False, faces_idx = 'all')
    - lines: bool
    - faces: bool
    - faces_idx: list of ints OR string
        Plots folded form. Plots edges if "lines" is true, plots faces if
        "faces" is true, and "faces_idx" are the indexes of the desired faces
        
    create_data_dict()
        Gets all of the variables on the origami and creates an up to date
        dict of all of it.
        
    save(outfile)
    - outfile: string
        calls create_data_dict() and then saves to path "outfile" in JSON 
        format
    
    Static Methods
    ---------
    load(infile)
    - infile: string
        loads FOLD file and creates an origami object with it


        """
    
    def __init__(self,frame_title = 'Origami fold file'):
        self.file_version = 1
        self.file_creator = 'FOLD_generator.py'
        self.file_author = 'Evandro Bernardes'
        self.file_classes = 'singleModel'
        self.frame_classes = 'foldedForm'
        self.frame_attributes = ['3D','orientable']
        self.frame_unit = 'unit'
        self.frame_title = frame_title
        self.idx = {}
        self.vertices_coords = []
        self.vertices_names = []
        self.vertices_size = 0
        self.faces_vertices = []
        self.faceOrders = []
        self.faces_size = 0
        self.edges_vertices = []
        self.edges_length = []
        self.edges_assignment = []
        self.edges_foldAngle = []
    
    
    def add_vertex(self,name,coords):
        """
            Adds "coords" to self.vertices_coords, "name" to self.vertices_names 
            and create a new entry in self.idx using "name" to point to this new
            vertex. Also updates self.vertices_size
            
            Parameters
            ----------
            name : string
                Name of vertex.
            coords : triplet
                Position of vertex
    
            Returns
            -------
            None.
            
                """
        if type(coords) == tuple or type(coords) == list:
            coords = np.array(coords)
        elif type(coords) != np.ndarray:
            raise ValueError
        if type(name) != str:
            raise ValueError
        
        self.vertices_names.append(name)
        self.idx[name] = self.vertices_size
        self.vertices_size = self.vertices_size + 1
        self.vertices_coords.append(coords)
    
    def connect(self,name_1,name_2,assignment,angle = 'null'):
        """
            Gets indices of vertices named "name_1" and "name_2" and create an
            edge between these two vertices in self.edges_vertices. Then,
            adds "assignment" to self.edges_assignments
            
            Parameters
            ----------
            name_1 : string
                Name of first vertex.
            name_2 : string
                Name of second vertex
            assignement : string
                Type of edge: 'M', 'V', 'U', 'F' or 'U'
    
            Returns
            -------
            None.
                """
        idx1 = self.idx[name_1]
        idx2 = self.idx[name_2]
        self.edges_vertices.append((idx1,idx2))
        self.edges_assignment.append(assignment)
        if angle != 'null':
            angle = abs(angle) if assignment == 'V' else -abs(angle)
        self.edges_foldAngle.append(angle)
        
        # calculate length
        coords1 = self.vertices_coords[idx1]
        coords2 = self.vertices_coords[idx2]
        self.edges_length.append(np.linalg.norm(coords1-coords2))
        
    
    def add_face(self,vertices_names):
        """
            Gets indices of vertices with names in "vertices_names" and create
            a face with these vertices, adding it to self.faces_vertices. Updates
            self.faces_size
            
            Parameters
            ----------
            vertices_names : list of strings
                List of names of vertices
    
            Returns
            -------
            None.
                """
        if type(vertices_names) != list:
            raise ValueError
        if type(vertices_names[0]) != str:
            raise ValueError
        
        vertex_idxs = []
        for vertex in vertices_names:
            vertex_idxs.append(self.idx[vertex])
        self.faces_vertices.append(vertex_idxs)
        
        self.faces_size = self.faces_size + 1
        
    def plot(self,lines=True, faces = False, faces_idx = 'all'):
        """
             Plots folded form. Plots edges if "lines" is true, plots faces if
            "faces" is true, and "faces_idx" are the indexes of the desired faces
    
            Parameters
            ----------
            lines : bool, optional
                Plot edges?. The default is True.
            faces : bool, optional
                Plot faces?. The default is False.
            faces_idx : list of ints or string, optional
                Indices of faces to plot. The default is 'all'.
    
            Returns
            -------
            None.
    
                """
        ax = plt.axes(projection="3d")
        
        if lines:
            for assignement,line in zip(self.edges_assignment,self.edges_vertices):
                v1,v2 = line
                if assignement == 'B':
                    color = 'k'
                elif assignement == 'V':
                    color = 'b'
                elif assignement == 'M':
                    color = 'r'
                elif assignement == 'C':
                    color = 'g'
                elif assignement == 'U':
                    color = 'm'
                elif assignement == 'F':
                    color = 'c'
                else:
                    color = 'g'
                
                if color != '':
                    x1,y1,z1 = self.vertices_coords[v1]
                    x2,y2,z2 = self.vertices_coords[v2]
                    ax.plot([x1,x2], [y1,y2], [z1,z2], color)
                
        if faces:
            if faces_idx == 'all':
                faces_vertices = self.faces_vertices
            else:
                if type(faces_idx) != list:
                    faces_idx = [faces_idx]
                faces_vertices = [self.faces_vertices[i] for i in faces_idx]
            for face in faces_vertices:
                x = []
                y = []
                z = []
                for p in face:
                    x_,y_,z_ = self.vertices_coords[p]
                    x.append(x_)
                    y.append(y_)
                    z.append(z_)
                verts = [list(zip(x,y,z))]
                ax.add_collection3d(Poly3DCollection(verts))
        
        plt.show()
        
    def create_data_dict(self):
        """
            Gets all of the variables on the origami and creates an up to date
            dict of all of it.
    
            Returns
            -------
            None.
    
            """
        vertices_coords = [list(v) for v in self.vertices_coords]
        # for v in vertices_coords:
        #     print(v)
        return collections.OrderedDict([
                ['file_version', self.file_version],
                ['file_creator', self.file_creator],
                ['file_author', self.file_author],
                ['file_classes', self.file_classes],
                ['frame_classes', self.frame_classes],
                ['frame_attributes', self.frame_attributes],
                ['frame_unit', self.frame_unit],
                ['frame_title', self.frame_title],
                ['vertices_names', self.vertices_names],
                ['vertices_coords', vertices_coords],
                ['faces_vertices', self.faces_vertices],
                ['faceOrders', self.faceOrders],
                ['edges_assignment', self.edges_assignment],
                ['edges_length', self.edges_length],
                # ['edges_foldAngle', self.edges_foldAngle],
                ['edges_vertices', self.edges_vertices]])
        
    def save(self,outfile):
        """
            Calls create_data_dict() and then saves to path "outfile" in JSON 
            format
    
            Parameters
            ----------
            outfile : string
                Path to output file.
    
            Returns
            -------
            None.
    
                """
        
        json_file = json.dumps(self.create_data_dict(),sort_keys=False,indent=2)
        f = open(outfile,'w')
        f.write(json_file)
        f.close()
        
    @staticmethod
    def load(infile):
        """
            Loads FOLD file and creates an origami object with it
    
            Parameters
            ----------
            infile : string
                Path to input file.
    
            Returns
            -------
            None.
    
                """
        with open(infile) as f:
            infile = json.load(f)
            
        origami = Origami()
        
        for key in infile.keys():
            command = "origami."+key+" = "+"infile['"+key+"']"
            # print(command)
            exec(command)
        
        size = len(origami.vertices_coords)
        origami.vertices_coords = [np.array(v) for v in origami.vertices_coords]
        
        # check possible problems 
        try:
            for i in range(len(origami.edges_vertices)):
                v1,v2 = origami.edges_vertices[i]
                if v1 > size or v2 > size:
                    print('Out of vertices range! I will try to fix it...')
                    v1 = v1 % origami.vertices_size
                    v2 = v2 % origami.vertices_size
                # print(str(v1)+' '+str(v2))
                origami.edges_vertices[i] = (int(v1),int(v2))
        except:
            pass
        
        # check possible problems 
        try:
            for i in range(len(origami.faces_vertices)):
                v1,v2,v3 = origami.faces_vertices[i]
                if v1 > size or v2 > size or v3 > size:
                    print('Out of vertices range! I will try to fix it...')
                    v1 = v1 % origami.vertices_size
                    v2 = v2 % origami.vertices_size
                    v3 = v3 % origami.vertices_size
                # print(str(v1)+' '+str(v2)+' '+str(v3))
                origami.faces_vertices[i] = [int(v1),int(v2),int(v3)]
                # print(origami.faces_vertices[i])
        except:
            pass
        
        origami.vertices_size = size
        
        # if edge lenghts are not found, calculate them
        if not hasattr(origami, 'edges_length'):
            origami.edges_length = []
            for idx1,idx2 in origami.edges_vertices:
                coords1 = origami.vertices_coords[idx1]
                coords2 = origami.vertices_coords[idx2]
                origami.edges_length.append(np.linalg.norm(coords1-coords2))
        
        return origami
    



    
