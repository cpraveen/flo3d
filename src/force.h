#ifndef __FORCE_H__
#define __FORCE_H__

struct ForceData
{
   std::string name;
   std::vector<int> face_type;
};

struct Force
{
   std::vector<int> face;
   Vector           value;
};

#endif
