/* \author Aaron Brown */
// Quiz on implementing kd tree

#include "../../render/render.h"


// Structure to represent node of kd tree
struct Node
{
	std::vector<float> point;
	int id;
	Node* left;
	Node* right;

	Node(std::vector<float> arr, int setId)
	:	point(arr), id(setId), left(NULL), right(NULL)
	{}
};

struct KdTree
{
	Node* root;

	KdTree()
	: root(NULL)
	{}

	void insertHelper(Node** node,int iteration,std::vector<float> point, int id){
		if(*node == NULL)
			*node = new Node(point,id);
		else
		{
			uint cd = iteration % point.size();

			if(point[cd] < (*node)->point[cd])
				insertHelper(&(*node)->left,++iteration, point, id);
			else
				insertHelper(&(*node)->right,++iteration, point, id);
			
		}
			
	}

	void insert(std::vector<float> point, int id)
	{
		// TODO: Fill in this function to insert a new point into the tree
		// the function should create a new node and place correctly with in the root
		insertHelper(&root,0,point,id);

	}

	void searchHelper(std::vector<float> target, Node* node, int iteration, float distanceTol, std::vector<int>& ids ){
		
		if(node != NULL){
			if(target.size() == 2){
				if((node->point[0] >= target[0]-distanceTol && node->point[0] <= target[0]+distanceTol) &&
			   (node->point[1] >= target[1]-distanceTol && node->point[1] <= target[1]+distanceTol)){
				

				float distance = sqrt(((node->point[0] - target[0])*(node->point[0] - target[0]))+
			   						 ((node->point[1] - target[1])*(node->point[1] - target[1])));
				if (distance <= distanceTol)
					ids.push_back(node->id);
			}
			}
			else{
				if((node->point[0] >= target[0]-distanceTol && node->point[0] <= target[0]+distanceTol) &&
			   (node->point[1] >= target[1]-distanceTol && node->point[1] <= target[1]+distanceTol) &&
			   (node->point[2] >= target[2]-distanceTol && node->point[2] <= target[2]+distanceTol)
			   ){
				

				float distance = sqrt(((node->point[0] - target[0])*(node->point[0] - target[0]))+
			   						 ((node->point[1] - target[1])*(node->point[1] - target[1])) +
			   						 ((node->point[2] - target[2])*(node->point[2] - target[2])));
				if (distance <= distanceTol)
					ids.push_back(node->id);
			}

			}
			

			if((target[iteration % target.size()]-distanceTol) <= node->point[iteration % target.size()])
				searchHelper(target, node->left, iteration + 1, distanceTol, ids );

			if((target[iteration % target.size()]+distanceTol) > node->point[iteration % target.size()])
				searchHelper(target, node->right, iteration + 1, distanceTol, ids );
			   
		}
		
	}

	// return a list of point ids in the tree that are within distance of target
	std::vector<int> search(std::vector<float> target, float distanceTol)
	{
		std::vector<int> ids;
		searchHelper(target, root, 0, distanceTol, ids );
		return ids;
	}
	

};




