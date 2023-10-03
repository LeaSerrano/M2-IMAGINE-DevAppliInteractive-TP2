#ifndef SKELETON_H
#define SKELETON_H


#include <vector>
#include <queue>
#include <map>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include "Vec3.h"

#include <cmath>

#include <GL/glut.h>


// -------------------------------------------
// Basic Skeleton class
// -------------------------------------------

struct Articulation {
    // membres :
    Vec3 p; // une position

    int fatherBone; // there should be only 1
    std::vector< unsigned int > childBones;

    void setFatherBone( int f ) {
        if( fatherBone >= 0 ) {
            assert(fatherBone == f);
        }
        fatherBone = f;
    }
    bool isRoot() const {
        return fatherBone == -1;
    }

    Articulation() : fatherBone(-1) { childBones.clear(); }
};

struct Bone {
    // membres :
    unsigned int joints[2];

    int fatherBone; // there should be only 1
    std::vector< unsigned int > childBones;

    void setFatherBone( int f ) {
        if( fatherBone >= 0 ) {
            assert(fatherBone == f);
        }
        fatherBone = f;
    }
    bool isRoot() const {
        return fatherBone == -1;
    }

    Bone() : fatherBone(-1) { childBones.clear(); }
};

struct BoneTransformation{
    // membres :
    Mat3 localRotation;

    Mat3 world_space_rotation;
    Vec3 world_space_translation;

    BoneTransformation() : localRotation(Mat3::Identity()) ,  world_space_rotation(Mat3::Identity()) , world_space_translation(0,0,0) {}
};

struct SkeletonTransformation{
    // membres :
    std::vector< BoneTransformation > bone_transformations;
    std::vector< Vec3 > articulations_transformed_position;

    void resize( unsigned int n_bones , unsigned int n_articulations ) {
        bone_transformations.resize( n_bones );
        articulations_transformed_position.resize( n_articulations );
    }
};


struct Skeleton {
    // membres :
    std::vector< Articulation > articulations;
    std::vector< Bone > bones;
    std::vector< unsigned int > ordered_bone_indices; // process them by order in the hierarchy

    void buildStructure() {
        ordered_bone_indices.clear();
        std::vector< unsigned int > rootBones; // why not have several

        for( unsigned int b = 0 ; b < bones.size() ; ++b ) {
            Articulation & a0 = articulations[ bones[b].joints[0] ];
            Articulation & a1 = articulations[ bones[b].joints[1] ];
            a0.childBones.push_back( b );
            a1.setFatherBone( b );
        }

        for( unsigned int aIdx = 0 ; aIdx < articulations.size() ; ++aIdx ) {
            Articulation & a = articulations[ aIdx ];
            if( a.isRoot() ) {
                for( unsigned int bIt = 0 ; bIt < a.childBones.size() ; ++bIt ) {
                    unsigned int b = a.childBones[bIt];
                    rootBones.push_back( b );
                }
            }
            else {
                unsigned int bfIdx = a.fatherBone;
                Bone & bf = bones[bfIdx];
                for( unsigned int bIt = 0 ; bIt < a.childBones.size() ; ++bIt ) {
                    unsigned int bcIdx = a.childBones[bIt];
                    Bone & bc = bones[bcIdx];
                    bc.setFatherBone( bfIdx );
                    bf.childBones.push_back( bcIdx );
                }
            }
        }

        for( unsigned int rIt = 0 ; rIt < rootBones.size() ; ++rIt ) {
            unsigned int rootboneIdx = rootBones[rIt];
            std::queue< unsigned int > bonesIndices;
            bonesIndices.push(rootboneIdx);
            while( ! bonesIndices.empty()) {
                unsigned int bIdx = bonesIndices.front();
                bonesIndices.pop();
                ordered_bone_indices.push_back(bIdx);
                Bone & b = bones[bIdx];
                for( unsigned int bIt = 0 ; bIt < b.childBones.size() ; ++bIt ) {
                    unsigned int bcIdx = b.childBones[bIt];
                    bonesIndices.push(bcIdx);
                }
            }
        }

        assert( ordered_bone_indices.size() == bones.size() );
    }

    void load (const std::string & filename){
        std::ifstream in (filename.c_str ());
        if (!in)
            exit (EXIT_FAILURE);
        std::string tmpString;
        unsigned int sizeA;
        in >> tmpString >> sizeA;
        articulations.resize (sizeA);
        for (unsigned int i = 0; i < sizeA; i++)
            in >> articulations[i].p[0] >> articulations[i].p[1] >> articulations[i].p[2];

        unsigned int sizeB;
        in >> tmpString >> sizeB;
        bones.resize (sizeB);
        for (unsigned int i = 0; i < sizeB; i++) {
            for (unsigned int j = 0; j < 2; j++)
                in >> bones[i].joints[j];
        }
        in.close ();

        buildStructure();
    }


    void computeGlobalTransformationParameters( SkeletonTransformation & transfo ) {
        std::vector< Vec3 > & articulations_transformed_position = transfo.articulations_transformed_position;
        articulations_transformed_position.resize( articulations.size() );
        for( unsigned int bIt = 0 ; bIt < ordered_bone_indices.size() ; ++bIt ) {
            unsigned bIdx = ordered_bone_indices[bIt];
            Bone & b = bones[bIdx];

            if( b.isRoot() ) {
                Vec3 a0RestPos = articulations[ b.joints[0] ].p;
                Vec3 a0TargetPos = a0RestPos;
                articulations_transformed_position[ b.joints[0] ] = a0TargetPos;
                BoneTransformation & bone_transformation = transfo.bone_transformations[bIdx];
                bone_transformation.world_space_rotation = bone_transformation.localRotation;

                // set the articulation as pivot point :
                bone_transformation.world_space_translation = a0TargetPos - bone_transformation.world_space_rotation * a0RestPos;

                // update the child articulation :
                Vec3 a1RestPos = articulations[ b.joints[1] ].p;
                Vec3 a1TargetPos = bone_transformation.world_space_rotation * a1RestPos + bone_transformation.world_space_translation;
                articulations_transformed_position[ b.joints[1] ] = a1TargetPos;
            }
            else{
                Vec3 a0RestPos = articulations[ b.joints[0] ].p;
                Vec3 a0TargetPos = articulations_transformed_position[ b.joints[0] ];

                BoneTransformation & bone_transformation = transfo.bone_transformations[bIdx];
                BoneTransformation & bFatherTransfo = transfo.bone_transformations[b.fatherBone];
                bone_transformation.world_space_rotation = bFatherTransfo.world_space_rotation * bone_transformation.localRotation;

                // set the articulation as pivot point :
                bone_transformation.world_space_translation = a0TargetPos - bone_transformation.world_space_rotation * a0RestPos;

                // update the child articulation :
                Vec3 a1RestPos = articulations[ b.joints[1] ].p;
                Vec3 a1TargetPos = bone_transformation.world_space_rotation * a1RestPos + bone_transformation.world_space_translation;
                articulations_transformed_position[ b.joints[1] ] = a1TargetPos;
            }
        }
    }




    void computeProceduralAnim( double t , SkeletonTransformation & transfo ) {
        transfo.bone_transformations.resize( bones.size() );
        for( unsigned int bIt = 0 ; bIt < ordered_bone_indices.size() ; ++bIt ) {
            unsigned bIdx = ordered_bone_indices[bIt];
            Bone & b = bones[bIdx];
            if( b.isRoot() ) {
                BoneTransformation & bone_transformation = transfo.bone_transformations[bIdx];
                bone_transformation.localRotation = Mat3::Identity();
            }
            else{
                BoneTransformation & bone_transformation = transfo.bone_transformations[bIdx];
                Vec3 axis( cos( 2 * M_PI * bIt / (double)(bones.size()) )  ,  sin( 2 * M_PI * bIt / (double)(bones.size()) )  , 0.0 );
                bone_transformation.localRotation = Mat3::getRotationMatrixFromAxisAndAngle( axis , (0.25*M_PI) * cos( t ) );
            }
        }

        // update articulation positions:
        computeGlobalTransformationParameters(transfo);
    }




//Version 1
/*void updateIKChain(SkeletonTransformation &transfoIK, unsigned int targetArticulationIdx, Vec3 targetPosition, unsigned int nombreIterations = 20, double epsilonPrecision = 0.000001) {
    for (unsigned int iteration = 0; iteration < nombreIterations; ++iteration) {
        BoneTransformation targetTransformation = transfoIK.bone_transformations[targetArticulationIdx];
        Articulation targetArtic = articulations[targetArticulationIdx];

        Vec3 error = transfoIK.articulations_transformed_position[targetArticulationIdx] - targetPosition;

        if (error.length() < epsilonPrecision) {
            return;
        }

        for (int i = ordered_bone_indices.size() - 1; i >= 0; --i) {
            Articulation currentArtic = articulations[i];

            Vec3 jointToTarget = targetPosition - currentArtic.p;

            Mat3 rotation = Mat3::getRotationMatrixAligning(jointToTarget, error);

            targetTransformation.localRotation = rotation * targetTransformation.localRotation;
            transfoIK.articulations_transformed_position[i] = currentArtic.p;

            error = targetPosition - currentArtic.p;

            if (error.length() < epsilonPrecision) {
                return;
            }
        }
    }
}*/


//Version 2
void updateIKChain(SkeletonTransformation &transfoIK, unsigned int targetArticulation, Vec3 targetPosition, unsigned int maxIterNumber = 20, double epsilonPrecision = 0.000001) {
    
    unsigned int numIterations = 0;

    while (numIterations < maxIterNumber) {
        for (int i = ordered_bone_indices.size() - 1; i >= 0; --i) {
            unsigned int boneIdx = ordered_bone_indices[i];
            Bone &bone = bones[boneIdx];
            unsigned int jointIdx0 = bone.joints[0];
            unsigned int jointIdx1 = bone.joints[1];
            Articulation &articulation0 = articulations[jointIdx0];
            Articulation &articulation1 = articulations[jointIdx1];

            Vec3 toTarget = targetPosition - articulation0.p;
            toTarget.normalize();

            Vec3 boneDirection = articulation1.p - articulation0.p;
            boneDirection.normalize();

            Mat3 rotationMatrix = Mat3::getRotationMatrixAligning(boneDirection, toTarget);

            BoneTransformation &boneTransfo = transfoIK.bone_transformations[boneIdx];
            boneTransfo.localRotation = rotationMatrix;

            for (int j = i; j < ordered_bone_indices.size(); ++j) {
                unsigned int updatedBoneIdx = ordered_bone_indices[j];
                Bone &updatedBone = bones[updatedBoneIdx];
                BoneTransformation &updatedBoneTransfo = transfoIK.bone_transformations[updatedBoneIdx];

                if (!updatedBone.isRoot()) {
                    unsigned int parentBoneIdx = updatedBone.fatherBone;
                    BoneTransformation &parentBoneTransfo = transfoIK.bone_transformations[parentBoneIdx];

                    updatedBoneTransfo.world_space_rotation = parentBoneTransfo.world_space_rotation * updatedBoneTransfo.localRotation;
                    updatedBoneTransfo.world_space_translation = parentBoneTransfo.world_space_translation + (parentBoneTransfo.world_space_rotation * boneDirection);
                } else {
                    updatedBoneTransfo.world_space_rotation = updatedBoneTransfo.localRotation;
                    updatedBoneTransfo.world_space_translation = articulation0.p;
                }
            }

            Articulation &endEffector = articulations[targetArticulation];
            Vec3 newEndEffectorPosition = endEffector.p;

            if ((targetPosition - newEndEffectorPosition).length() < epsilonPrecision) {
                return; 
            }
        }

        numIterations++;
    }
}

    //----------------------------------------------//
    //----------------------------------------------//
    //----------------------------------------------//
    // draw functions :
    //----------------------------------------------//
    //----------------------------------------------//
    //----------------------------------------------//
    void draw( int displayedBone , int displayedArticulation ) const {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);
        glLineWidth(3.0);
        glBegin (GL_LINES);
        for (unsigned int i = 0; i < bones.size (); i++) {
            glColor3f(1,0,0);
            {
                const Articulation & v = articulations[bones[i].joints[0]];
                glVertex3f (v.p[0], v.p[1], v.p[2]);
            }
            glColor3f(1,1,1);
            {
                const Articulation & v = articulations[bones[i].joints[1]];
                glVertex3f (v.p[0], v.p[1], v.p[2]);
            }
        }
        glEnd ();

        // we highlight the ordered bone number displayedBone
        if( displayedBone >= 0 && displayedBone < ordered_bone_indices.size() ) {
            displayedBone = ordered_bone_indices[displayedBone];
            glLineWidth(8.0);
            glBegin (GL_LINES);
            glColor3f(1,0,0);
            {
                const Articulation & v = articulations[bones[displayedBone].joints[0]];
                glVertex3f (v.p[0], v.p[1], v.p[2]);
            }
            glColor3f(1,0,0);
            {
                const Articulation & v = articulations[bones[displayedBone].joints[1]];
                glVertex3f (v.p[0], v.p[1], v.p[2]);
            }
            glEnd ();
        }

        // draw articulations:
        glPointSize(12.0);
        glBegin(GL_POINTS);
        glColor3f(0.5,0,0);
        for (unsigned int i = 0; i < articulations.size (); i++) {
            const Articulation & v = articulations[i];
            glVertex3f (v.p[0], v.p[1], v.p[2]);
        }
        glEnd();

        if( displayedArticulation >= 0 && displayedArticulation < articulations.size() ) {
            glPointSize(16.0);
            glBegin(GL_POINTS);
            glColor3f(1,0,0);
            {
                const Articulation & v = articulations[displayedArticulation];
                glVertex3f (v.p[0], v.p[1], v.p[2]);
            }
            glEnd();
        }

        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);
    }
    void drawTransformedSkeleton( int displayedBone , int displayedArticulation , SkeletonTransformation const & transfo ) const {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);
        glLineWidth(3.0);
        glBegin (GL_LINES);
        for (unsigned int i = 0; i < bones.size (); i++) {
            glColor3f(1,0,0);
            {
                Vec3 p = transfo.articulations_transformed_position[ bones[i].joints[0] ];
                glVertex3f (p[0], p[1], p[2]);
            }
            glColor3f(1,1,1);
            {
                Vec3 p = transfo.articulations_transformed_position[ bones[i].joints[1] ];
                glVertex3f (p[0], p[1], p[2]);
            }
        }
        glEnd ();

        // we highlight the ordered bone number displayedBone
        if( displayedBone >= 0 && displayedBone < ordered_bone_indices.size() ) {
            displayedBone = ordered_bone_indices[displayedBone];
            glLineWidth(8.0);
            glBegin (GL_LINES);
            glColor3f(1,0,0);
            {
                Vec3 p = transfo.articulations_transformed_position[ bones[displayedBone].joints[0] ];
                glVertex3f (p[0], p[1], p[2]);
            }
            glColor3f(1,0,0);
            {
                Vec3 p = transfo.articulations_transformed_position[ bones[displayedBone].joints[1] ];
                glVertex3f (p[0], p[1], p[2]);
            }
            glEnd ();
        }

        // draw articulations:
        glPointSize(12.0);
        glBegin(GL_POINTS);
        glColor3f(0.5,0,0);
        for (unsigned int i = 0; i < articulations.size (); i++) {
            Vec3 p = transfo.articulations_transformed_position[ i ];
            glVertex3f (p[0], p[1], p[2]);
        }
        glEnd();

        if( displayedArticulation >= 0 && displayedArticulation < articulations.size() ) {
            glPointSize(16.0);
            glBegin(GL_POINTS);
            glColor3f(1,0,0);
            {
                Vec3 p = transfo.articulations_transformed_position[ displayedArticulation ];
                glVertex3f (p[0], p[1], p[2]);
            }
            glEnd();
        }

        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);
    }
};



#endif // SKELETON_H
