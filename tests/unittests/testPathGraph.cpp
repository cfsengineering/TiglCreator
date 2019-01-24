//
// Created by cfse on 1/24/19.
//

#include "test.h"

#include "CTiglError.h"
#include "PathGraph.h"
#include "tigl.h"

#include "CCPACSConfiguration.h"
#include "CCPACSConfigurationManager.h"
#include "CCPACSFuselage.h"
#include "CCPACSFuselages.h"

#include <string.h>

TEST(TestPathGraph, emptyGraph) {

  tigl::PathGraph emptyGraph;

  EXPECT_TRUE(emptyGraph.isEmpty());
  EXPECT_TRUE(emptyGraph.getGraphAsVector().size() == 0);
  EXPECT_TRUE(emptyGraph.getElementsInBetween("fsadf", "fdsaf").size() == 0);
  EXPECT_TRUE(emptyGraph.getElementsBefore("aa").size() == 0);
  EXPECT_TRUE(emptyGraph.getElementsAfter("fsa").size() == 0);
  EXPECT_FALSE(emptyGraph.contains("fsda"));
  EXPECT_THROW(emptyGraph.getRoot(), tigl::CTiglError);
  EXPECT_THROW(emptyGraph.getLeaf(), tigl::CTiglError);
};

TEST(TestPathGraph, simpleGraph) {

  std::vector<std::string> simpleVectorGraph, expectedR;

  simpleVectorGraph.push_back("aa");
  simpleVectorGraph.push_back("bb");
  simpleVectorGraph.push_back("cc");
  simpleVectorGraph.push_back("dd");
  simpleVectorGraph.push_back("ee");
  simpleVectorGraph.push_back("ff");
  simpleVectorGraph.push_back("gg");

  tigl::PathGraph simpleGraph;

  simpleGraph.set(simpleVectorGraph);

  EXPECT_FALSE(simpleGraph.isEmpty());

  EXPECT_TRUE(simpleGraph.getGraphAsVector() == simpleVectorGraph);

  expectedR.clear();
  expectedR.push_back("cc");
  expectedR.push_back("dd");
  expectedR.push_back("ee");
  EXPECT_TRUE(simpleGraph.getElementsInBetween("cc", "ee") == expectedR);
  EXPECT_TRUE(simpleGraph.getElementsInBetween("ee", "cc") == expectedR);

  EXPECT_TRUE(simpleGraph.getElementsInBetween("fsadf", "ee").size() == 0);

  expectedR.clear();
  expectedR.push_back("aa");
  expectedR.push_back("bb");
  EXPECT_TRUE(simpleGraph.getElementsBefore("cc") == expectedR);

  EXPECT_TRUE(simpleGraph.getElementsBefore("aa").size() == 0);

  expectedR.clear();
  expectedR.push_back("ff");
  expectedR.push_back("gg");
  EXPECT_TRUE(simpleGraph.getElementsAfter("ee") == expectedR);

  EXPECT_TRUE(simpleGraph.getElementsAfter("fs").size() == 0);

  EXPECT_FALSE(simpleGraph.contains("fsda"));
  EXPECT_TRUE(simpleGraph.contains("dd"));

  EXPECT_TRUE(simpleGraph.getRoot() == "aa");
  EXPECT_TRUE(simpleGraph.getLeaf() == "gg");
}

TEST(TestPathGraph, malformedGraph) {

  std::vector<std::string> malformedVectorGraph;

  malformedVectorGraph.push_back("aa");
  malformedVectorGraph.push_back("bb");
  malformedVectorGraph.push_back("bb");

  tigl::PathGraph graph;

  EXPECT_THROW(graph.set(malformedVectorGraph), tigl::CTiglError);
}

TEST(TestFuselageGraph, simpleGraph) {

  // open cpacs file

  TiglCPACSConfigurationHandle tiglHandle = -1;
  TixiDocumentHandle tixiHandle = -1;

  std::string filename = "TestData/simpletest.cpacs.xml";
  ASSERT_EQ(SUCCESS, tixiOpenDocument(filename.c_str(), &tixiHandle));
  ASSERT_EQ(TIGL_SUCCESS,
            tiglOpenCPACSConfiguration(tixiHandle, "", &tiglHandle));

  tigl::CCPACSConfigurationManager *manager = nullptr;
  tigl::CCPACSConfiguration *config = nullptr;
  tigl::CCPACSFuselage *fuselage = nullptr;

  manager = &(tigl::CCPACSConfigurationManager::GetInstance());
  config = &(manager->GetConfiguration(tiglHandle));
  fuselage = &(config->GetFuselage(1));

  tigl::FuselageGraph fuselageGraph;
  fuselageGraph.set(fuselage);

  std::vector<std::string> expectedR;
  expectedR.push_back("D150_Fuselage_1Section1IDElement1");
  expectedR.push_back("D150_Fuselage_1Section2IDElement1");
  expectedR.push_back("D150_Fuselage_1Section3IDElement1");

  EXPECT_TRUE(fuselageGraph.getGraphAsVector() == expectedR);
}

TEST(TestFuselageGraph, invalidInput) {
  tigl::FuselageGraph fuselageGraph;
  EXPECT_THROW(fuselageGraph.set(nullptr), tigl::CTiglError);
}
