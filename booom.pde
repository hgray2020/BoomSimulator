
ArrayList<Mover> movers = new ArrayList<Mover>();
PVector gravity = new PVector(0, 0.2);
PVector groundBounce = new PVector(0, 1000);
int groundY = 0;
float explodePow = 10027000;
float friction = 0.1;
float bounceC = -0;
float bMass = 15;
float timescale = 0.01; 

void setup() {
  size(1500, 1500); 
  background(255, 0, 0);
  for (int i = 0; i < 75; i++) {
    movers.add(new Mover(new PVector(((i% (width/((2*bMass)+10)))*((2*bMass)+10)), 1200+((i/10)*((2*bMass)+10))), 200, false));
  }

  groundY = height-100;
}

void mouseClicked() {
  for (Mover m : movers) {
    float a = atan2(m.location.y-mouseY, m.location.x-mouseX);
    float d = dist(m.location.x, m.location.y, mouseX, mouseY);
    float forcex = (explodePow*cos(a))/d;
    float forcey = (explodePow*sin(a))/d;
    m.applyForce(new PVector(forcex, forcey));
  }
  //movers.add(new Mover(new PVector(750, -500), 500, true));
}
void keyReleased() {
}

void draw() {
  background(255); 
  for (int p = 0; p < movers.size(); p++) {
    Mover m = movers.get(p);
    m.update(); 
    m.display();
  }
  for (int p = 0; p < movers.size(); p++) {
    Mover m = movers.get(p);
    if (m.broke) {
      movers.remove(p);
    }
  }
  for (Mover m : movers) {
    if (m.ast) {
      for (Mover m2 : movers) {
        if (m.location.x != m2.location.x && m.location.y != m2.location.y) {
          //m.checkcollision(m2);
        }
      }
    }
  }
  if (keyPressed && key == 'w') {
    for (Mover m : movers) {
      m.velocity.mult(0);
    }
  }
  if (keyPressed) {
    if (key == 'q') {
      timescale = 0.01;
    }
    if (key == 'e') {
      timescale = 1;
    }
    if (key == 't') {
      timescale = 2;
    }
    if (key == 'y') {
      timescale = 4;
    }
  }



  /*
  for (int i = 0; i < movers.size(); i++) {
   Mover m = movers.get(i);
   for (int j = i+1; j < movers.size(); j++) {
   if (i != j) {
   m.checkCollision(movers.get(j));
   }
   }
   }*/
  //println(movers.get(0).velocity.x);
  for (int i = 0; i < movers.size(); i++) {
    Mover m = movers.get(i);
    for (int j = i+1; j < movers.size(); j++) {
      m.checkcollision(movers.get(j));
    }
  }
}



class Mover {
  PVector location; 
  PVector velocity; 
  PVector acceleration; 
  float mass; 
  int col = color(0, 0, 255);
  boolean ast = false;
  boolean broke = false;

  Mover(PVector l, float m, boolean asteroid) {
    location = l; 
    mass = m; 
    velocity = new PVector(0, 0); 
    acceleration = new PVector(0, 0);
    ast = asteroid;
    if (ast) {
      col = color(130);
    }
  }


  void update() {
    if (location.y < groundY-(mass/2) || dist(location.x, location.y, 750, groundY) < 250+(mass/2)) {

      acceleration.add(gravity);
    } else {
      if (velocity.x > friction) {
        velocity.x -= friction;
      } else if (velocity.x < -1*friction) {
        velocity.x += friction;
      } else {
        velocity.x = 0;
      }
    }

    bounce(); 
    velocity.add(acceleration); 
    if (acceleration.mag() > mass/8 && mass > 50) {
      float num = 20;
      float aHeading = velocity.heading();
      float c = sqrt(8);
      for (int i = 0; i < num; i++) {
        float a1;
        if (i < num/2) {
          a1 = (aHeading-(PI/6))+((PI/(num/2))*i);
        } else {
          a1 = aHeading+(PI/6+((PI/(num/2))*(i-(num/2))));
        }
        if (a1 < -PI) {
          a1+=TWO_PI;
        }
        if (a1 > PI) {
          a1-=TWO_PI;
        }
        Mover m;
        float v = this.velocity.mag()*c;
        if (i < num/2) {
          m = new Mover(new PVector(location.x-((mass/2)*cos(aHeading)), location.y-((mass/2)*sin(aHeading))), mass/num, false); //(random(num-((3*num)/8), num+((3*num)/8)))
        } else {
          m = new Mover(new PVector(location.x-((mass/2)*cos(aHeading)), location.y-((mass/2)*sin(aHeading))), mass/num, false);
        }
        
        
        PVector v1 = new PVector(v*cos(a1), v*sin(a1));
        
        m.velocity = v1;
        
        movers.add(m);
      } 
      broke = true;
    }

    location.add(new PVector(velocity.x*timescale, velocity.y*timescale)); 
    acceleration.mult(0);
  }
  void display() {
    noStroke();
    fill(col); 
    //rect(location.x-(mass/2), location.y-(mass/2), mass, mass);
    ellipse(location.x, location.y, mass, mass);
  }
  void applyForce(PVector force) {

    force.div(mass); 

    acceleration.add(force);
    println(acceleration.mag());
  }

  void bounce() {
    if (location.y > groundY-(mass/2)) {

      location.y = groundY-(mass/2); 
      velocity.y = -0.6*velocity.y; 

      // applyForce(groundBounce);
    }
    if (location.y < (mass/2)) {

      //location.y = (mass/2); 
      //velocity.y = -0.6*velocity.y; 

      // applyForce(groundBounce);
    }
    if (location.x > width-(mass/2)) {

      location.x = width-(mass/2); 
      velocity.x = -0.6*velocity.x;


      applyForce(new PVector(0.3, 0));
    }
    if (location.x < (mass/2)) {

      location.x = (mass/2); 
      velocity.x = -0.6*velocity.x;


      applyForce(new PVector(0.3, 0));
    }
  }
  void checkcollision(Mover other) {
    float ang = atan2(other.location.y-location.y, other.location.x-location.x);
    float d = dist(location.x, location.y, other.location.x, other.location.y);
    float mr = mass/(mass+other.mass);
    float omr = other.mass/(mass+other.mass);
    if (ast || other.ast) {
      if (other.ast) {
        println(omr);
      }
    }
    if (d <= (mass/2)+(other.mass/2)) {
      if (!other.ast) {

        float pow = 2;

        PVector oloc = new PVector(location.x+(mass*cos(ang)), location.y+(mass*sin(ang)));
        PVector loc = new PVector(other.location.x-(other.mass*cos(ang)), other.location.y-(other.mass*sin(ang)));
        location = loc;
        if (!ast) {
          other.location = oloc;
        }

        other.velocity.x = (mr*cos(ang)*mag(other.velocity.x, other.velocity.y))+(omr*cos(ang)*mag(velocity.x, velocity.y));
        other.velocity.y = (mr*sin(ang)*mag(other.velocity.x, other.velocity.y))+(omr*sin(ang)*mag(velocity.x, velocity.y));
        velocity.x = -1*((omr*cos(ang)*mag(velocity.x, velocity.y))+(mr*cos(ang)*mag(other.velocity.x, other.velocity.y)));
        velocity.y = -1*((omr*sin(ang)*mag(velocity.x, velocity.y))+(mr*sin(ang)*mag(other.velocity.x, other.velocity.y)));

        other.velocity.mult(1);
        velocity.mult(1);

        /*if (d < mass*0.5) {
         other.velocity.x = pow*(mr*cos(ang)*mag(other.velocity.x, other.velocity.y))+(omr*cos(ang)*mag(velocity.x, velocity.y));
         other.velocity.y = pow*(mr*sin(ang)*mag(other.velocity.x, other.velocity.y))+(omr*sin(ang)*mag(velocity.x, velocity.y));
         velocity.x = -1*pow*((omr*cos(ang)*mag(velocity.x, velocity.y))+(mr*cos(ang)*mag(other.velocity.x, other.velocity.y)));
         velocity.y = -1*pow*((omr*sin(ang)*mag(velocity.x, velocity.y))+(mr*sin(ang)*mag(other.velocity.x, other.velocity.y)));
         }*/
      } else {
        if (other.ast) {
          velocity.x = -1*(omr*cos(ang)*mag(other.velocity.x, other.velocity.y))+(mr*cos(ang)*mag(velocity.x, velocity.y));
          velocity.y = -1*(omr*sin(ang)*mag(other.velocity.x, other.velocity.y))+(mr*sin(ang)*mag(velocity.x, velocity.y));

          PVector oloc = new PVector(location.x+(other.mass*cos(ang)), location.y+(mass*sin(ang)));
          PVector loc = new PVector(oloc.x-(other.mass*cos(ang)), oloc.y-(other.mass*sin(ang)));
          location = loc;
        }
      }
    }
  }
  void checkCollision(Mover other) {

    // Get distances between the balls components
    PVector distanceVect = PVector.sub(other.location, location);

    // Calculate magnitude of the vector separating the balls
    float distanceVectMag = distanceVect.mag();

    // Minimum distance before they are touching
    float minDistance = (mass + other.mass)/2;

    if (dist(location.x, location.y, other.location.x, other.location.y) < minDistance) {
      float distanceCorrection = (minDistance-distanceVectMag)/2.0;
      PVector d = distanceVect.copy();
      PVector correctionVector = d.normalize().mult(distanceCorrection);
      other.location.add(correctionVector);
      location.sub(correctionVector);

      // get angle of distanceVect
      float theta = distanceVect.heading();
      // precalculate trig values
      float sine = sin(theta);
      float cosine = cos(theta);

      /* bTemp will hold rotated ball locations. You 
       just need to worry about bTemp[1] location*/
      PVector[] bTemp = {
        new PVector(), new PVector()
      };

      /* this ball's location is relative to the other
       so you can use the vector between them (bVect) as the 
       reference point in the rotation expressions.
       bTemp[0].location.x and bTemp[0].location.y will initialize
       automatically to 0.0, which is what you want
       since b[1] will rotate around b[0] */
      bTemp[1].x  = cosine * distanceVect.x + sine * distanceVect.y;
      bTemp[1].y  = cosine * distanceVect.y - sine * distanceVect.x;

      // rotate Temporary velocities
      PVector[] vTemp = {
        new PVector(), new PVector()
      };

      vTemp[0].x  = cosine * velocity.x + sine * velocity.y;
      vTemp[0].y  = cosine * velocity.y - sine * velocity.x;
      vTemp[1].x  = cosine * other.velocity.x + sine * other.velocity.y;
      vTemp[1].y  = cosine * other.velocity.y - sine * other.velocity.x;

      /* Now that velocities are rotated, you can use 1D
       conservation of momentum equations to calculate 
       the final velocity along the x-axis. */
      PVector[] vFinal = {  
        new PVector(), new PVector()
      };

      // final rotated velocity for b[0]
      vFinal[0].x = ((mass - other.mass) * vTemp[0].x + 2 * other.mass * vTemp[1].x) / (mass + other.mass);
      vFinal[0].y = vTemp[0].y;

      // final rotated velocity for b[0]
      vFinal[1].x = ((other.mass - mass) * vTemp[1].x + 2 * mass * vTemp[0].x) / (mass + other.mass);
      vFinal[1].y = vTemp[1].y;

      // hack to avoid clumping
      bTemp[0].x += vFinal[0].x;
      bTemp[1].x += vFinal[1].x;

      /* Rotate ball locations and velocities back
       Reverse signs in trig expressions to rotate 
       in the opposite direction */
      // rotate balls
      PVector[] bFinal = { 
        new PVector(), new PVector()
      };

      bFinal[0].x = cosine * bTemp[0].x - sine * bTemp[0].y;
      bFinal[0].y = cosine * bTemp[0].y + sine * bTemp[0].x;
      bFinal[1].x = cosine * bTemp[1].x - sine * bTemp[1].y;
      bFinal[1].y = cosine * bTemp[1].y + sine * bTemp[1].x;

      // update balls to screen location
      other.location.x = location.x + bFinal[1].x;
      other.location.y = location.y + bFinal[1].y;

      location.add(bFinal[0]);

      // update velocities
      velocity.x = cosine * vFinal[0].x - sine * vFinal[0].y;
      velocity.y = cosine * vFinal[0].y + sine * vFinal[0].x;
      other.velocity.x = cosine * vFinal[1].x - sine * vFinal[1].y;
      other.velocity.y = cosine * vFinal[1].y + sine * vFinal[1].x;
    }
  }
}