#include <gtest/gtest.h>
#include "../../../src/Model/PhysicalProperties.hpp"
using namespace GlobalFlow::Model;

TEST(PhysicalProperty,AddAndSet){
	PropertyRepository<PhysicalProperty<int>> prop;
	prop.set<int>(5);
	ASSERT_EQ(prop.get<int>(),5);
	prop.set<int>(6);
	ASSERT_NE(prop.get<int>(),5);
}
TEST(PhysicalProperty,Emplace){
	struct test1;
	struct test2;
	PropertyRepository<PhysicalProperty<int,test1>,PhysicalProperty<int,test2>> prop;
	prop.emplace<int, test1>(5);
	prop.emplace<int,test2>(6);
	ASSERT_EQ((prop.get<int, test1>()), 5);
	ASSERT_EQ((prop.get<int, test2>()), 6);
	ASSERT_EQ((prop.get<int, test1>()), 5);
}
TEST(PhysicalProperty,AddTo){
	struct test1;
	struct test2;
	PropertyRepository<PhysicalProperty<int,test1>,PhysicalProperty<int,test2>> prop;
	ASSERT_THROW((prop.addTo<int, test1>(1)), std::invalid_argument);
	prop.set<int,test1>(5);
	prop.addTo<int,test1>(1);
	ASSERT_EQ((prop.get<int, test1>()), 6);
}
