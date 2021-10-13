use macroquad::prelude::*;
use rayon::prelude::*;
use roots;

const DT_PER_DISPLAY: f32 = 0.1; // how much time must pass (forward or backward) before displaying stuff
const DT_INIT: f32 = 0.1; // starting time-step
const C_INIT: f32 = 1.0; // e
const MAX_SLOWDOWN: f32 = 0.01; //0.01;// f32::MIN_POSITIVE; // dt*MAX_SLOWDOWN = minimum acceptable time-step
const FPS_WINDOW: usize = 50;

const GRAVITY_STRENGTH: f32 = 5.0;
const R_INIT: f32 = 10.0;
const INTEGRATION_METHOD_INIT: IntegrationMethod = IntegrationMethod::Euler;

const PI: f32 = std::f32::consts::PI;


#[derive(Copy, Clone)]
struct Disc{
    pos: Vec2,
    vel: Vec2,
    force: Vec2,
    r: f32,
    r2: f32,
    charge: f32,
    mass: f32,
    invmass: f32,
    id: i32
}

enum IntegrationMethod{
    Euler,
    Verlet,
}

impl Disc {
    fn new(x: f32, y: f32, vx: f32, vy: f32, r: f32, mass_density: f32, charge_density: f32, id: i32) -> Disc{
    Disc {pos: vec2(x,y),
          vel: vec2(vx, vy),
          force: vec2(0.0, 0.0),
          r: r,
          r2: r*r,
          charge: charge_density*2.0*PI*r*r, 
          mass: mass_density*2.0*PI*r*r,
          invmass: (1.0/(mass_density*2.0*PI*r*r)).abs(),
          id: id}
    }

    fn verlet_pos_step(&mut self, dt: &f32){
        self.pos += self.vel*(*dt) + 0.5*(self.force*self.invmass)*(*dt*dt);
    }

    fn verlet_vel_step(&mut self, new_force: &Vec2, dt: &f32){
        self.vel += 0.5*(self.force + *new_force)*self.invmass*(*dt);
        self.force = *new_force;
    }

    fn euler_pos_step(&mut self, dt: &f32){
        self.pos += self.vel*(*dt);
    }

    fn euler_vel_step(&mut self, new_force: &Vec2, dt: &f32){
        self.vel += self.force*self.invmass*(*dt);
        self.force = *new_force
    }

    fn time_until_collision_euler(&self, disc: &Disc) -> Option<f32>{
        let dx = self.pos - disc.pos;
        let dv = self.vel - disc.vel;

        // apply quadratic formula to solve for dt
        // || (x1-x2) + (v1-v2) || = r1 + r2
        // a2*dt^2 + a1*dt + a0 = 0

        // if dX dot dV is negative, entities traveling away from each other
        // so there cannot be a collision
        // (dt can only be positive in this situation if one disc already overlaps the other)

        let dv2 = dv.length_squared(); // a2, always positive
        let dx_dot_dv = dx.dot(dv); // equivalent to a1/2 in above Eq.
        let dx_sub_r2 = dx.length_squared() - (self.r + disc.r)*(self.r + disc.r); // a0

        // if discs don't intersect and are traveling away from each other,
        // there are no future collisions 
        // alternatively decarte's rule of signs implies no positive roots 
        // bc a2, a1, and a0 are the same sign
        if dx_sub_r2 >= 0.0 && dx_dot_dv >= 0.0 {
            None
        } else {
            // calculate 0.25*discriminant
            // discriminant = a1^2 - 4*a2*a0 = 4*( (1/4)*a1^2 - a2*a0 )
            let quarter_discriminant = dx_dot_dv * dx_dot_dv - dv2 * dx_sub_r2;

            if quarter_discriminant < 0.0 {
                // if solutions are imaginary, no collisions occur
                None
            } else {
                // calculate smallest positive solution
                // dt = ( -a1 - 2*sqrt(quarter_discriminant) )/ 2*a2
                //    = ( -2*dXdotdV - 2*sqrt(quarter_discriminant) )/ 2*a2
                //    = (-dXdotdV - sqrt(quarter_discriminant))/dV2 (because the twos cancel)
                let dt = (-dx_dot_dv - (quarter_discriminant).sqrt())/dv2;
                // if dt > 0.0 {Some(dt)} else {Some(0.0)}

                // if objects are traveling toward each other, its possible
                // for time until collision to be negative if the objects have
                // already penetrated. Handling these cases is crucial to prevent
                // free movement
                // assert!( (-dx_dot_dv + (quarter_discriminant).sqrt()) >= 0.0, "{0}, {1}", dt, (-dx_dot_dv + (quarter_discriminant).sqrt()));
                Some(dt)
            }
        }
    }

    /*
    fn time_until_collision_verlet(&self, disc: &Disc) -> Option<f32>{
        let dx = self.pos - disc.pos;
        let dv = self.vel - disc.vel;
        let da = self.force*self.invmass - disc.force*disc.invmass;

        // apply quartic formula to solve for dt
        let a4 = 0.25*da.length_squared(); // always positive
        let a3 = da.dot(dv);
        let a2 = dv.length_squared() + da.dot(dx);
        let a1 = 2.0*dv.dot(dx);
        let a0 = dx.length_squared() - (self.r + disc.r)*(self.r + disc.r);

        if a0 >= 0.0 && a1 >= 0.0 && a3 >= 0.0 {
            // if non-intersecting, velocity pointed away, acceleration pointed away
            // there will not be a collision
            None
        } else {
            let root_list = roots::find_roots_quartic(a4, a3, a2, a1, a0);
            if a0 >= 0.0 {
                // return smallest positive time
                match root_list {
                    roots::Roots::No(_) => None,
                    roots::Roots::One(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt > 0.0).into_iter().reduce(|a, b| if a<b {a} else {b})},
                    roots::Roots::Two(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt > 0.0).into_iter().reduce(|a, b| if a<b {a} else {b})},
                    roots::Roots::Three(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt > 0.0).into_iter().reduce(|a, b| if a<b {a} else {b})},
                    roots::Roots::Four(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt > 0.0).into_iter().reduce(|a, b| if a<b {a} else {b})}
                }
            } else {
                // return smallest-in-magnitude negative time
                match root_list {
                    roots::Roots::No(_) => None,
                    roots::Roots::One(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt < 0.0).into_iter().reduce(|a, b| if a>b {a} else {b})},
                    roots::Roots::Two(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt < 0.0).into_iter().reduce(|a, b| if a>b {a} else {b})},
                    roots::Roots::Three(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt < 0.0).into_iter().reduce(|a, b| if a>b {a} else {b})},
                    roots::Roots::Four(dts) => {std::array::IntoIter::new(dts).filter(|&dt| dt < 0.0).into_iter().reduce(|a, b| if a>b {a} else {b})}
                }
            }
        }
    }
    */

    fn time_until_collision(&self, disc: &Disc, method: &IntegrationMethod) -> Option<f32>{
        match method {
            IntegrationMethod::Euler => self.time_until_collision_euler(disc),
            IntegrationMethod::Verlet => self.time_until_collision_euler(disc)
        }
    }

    fn dist_vec(&self, disc: &Disc) -> Vec2{
        disc.pos - self.pos
    }

    fn distance(&self, disc: &Disc) -> f32{
        self.dist_vec(disc).length()
    }

    fn distance_recip(&self, disc: &Disc) -> f32{
        self.dist_vec(disc).length_recip()
    }

    fn draw(&self){
        draw_circle(self.pos.x, self.pos.y, self.r, BLACK);
        draw_circle(self.pos.x, self.pos.y, if self.r-1.0 > 0.0 {self.r-1.0} else {0.0}, WHITE);
        
        if self.mass > 0.0 {
            draw_text("+", self.pos.x-5.0, self.pos.y+5.0, 25.0, BLACK);
        }else if self.mass < -0.0 {
            draw_text("-", self.pos.x-5.0, self.pos.y+5.0, 25.0, BLACK);
        } else {
            draw_text("o", self.pos.x-5.0, self.pos.y+5.0, 25.0, BLACK);
        }
    }
}

fn get_random_unit_vec() -> Vec2{
    let x: f32 = rand::RandomRange::gen_range(0.0,1.0);
    let y = (1.0-x*x).sqrt();
    vec2(x, y)
}

fn inelastic_collision(disc1: Disc, disc2: Disc, c: f32) -> (Vec2, Vec2) {
    let dx = disc1.pos - disc2.pos;
    let dx_mag2 = dx.length_squared();
    if dx_mag2 > 0.0 {
        let dv = disc1.vel - disc2.vel;
        let m1 = disc1.mass.abs();
        let m2 = disc2.mass.abs();

        // calculate momentum transfer if there's a normal vector
        let m1m2 = m1+m2;
        let vdotx = dx.dot(dv);

        let term1 = (1.0+c)*(vdotx/dx_mag2)*dx/(m1+m2);
        (-m2*term1, m1*term1)
    } else {
        // otherwise don't perform momentum transfer
        (vec2(0.0, 0.0), vec2(0.0, 0.0))
    }
}

fn gravity(disc1: &Disc, disc2: &Disc) -> Vec2{
    let vec_diff = disc1.dist_vec(disc2);
    let distance = vec_diff.length();
    let d_avg = disc1.r + disc2.r;

    if distance > d_avg {
        GRAVITY_STRENGTH*disc1.mass*disc2.mass*vec_diff/(distance*distance*distance)
    } else if distance > 0.0 {
        - GRAVITY_STRENGTH*disc1.mass.abs()*disc2.mass.abs()*vec_diff/(d_avg*d_avg*distance)
    } else {
        - GRAVITY_STRENGTH*disc1.mass.abs()*disc2.mass.abs()*get_random_unit_vec()/(d_avg*d_avg)
    }       
}

fn force_sum(force: &dyn Fn(&Disc, &Disc) -> Vec2, disc1: &Disc, disc_list: &Vec<Disc>) -> Vec2{
    // 14 FPS, 400 entities, dt=0.001
    // disc_list.iter().filter(|disc2| disc1.id != disc2.id).map(|disc2| force(disc1, disc2)).collect::<Vec<Vec2>>().iter().sum()

    // 14 FPS, 400 entities, dt=0.001
    // disc_list.par_iter().filter(|disc2| disc1.id != disc2.id).map(|disc2| gravity_repulsion(disc1, disc2)).reduce(||vec2(0.0,0.0), |a,b| a+b)f

    // 30 FPS, 400 entities, dt=0.001
    disc_list.iter().filter(|disc2| disc1.id != disc2.id).map(|disc2| force(disc1, disc2)).fold(vec2(0.0,0.0), |a,b| a+b)
}

// Increase radius if W is pressed, decrease if S is pressed
fn update_r(r: &mut f32){
    if macroquad::input::is_key_down(macroquad::input::KeyCode::W){
        *r+=1.0;
    }
    if macroquad::input::is_key_down(macroquad::input::KeyCode::S){
        *r-=1.0;
    }
}

// Increase dt when D is pressed, decrease if A is pressed
fn update_dt(dt: &mut f32){
    let shift = macroquad::input::is_key_down(macroquad::input::KeyCode::LeftShift) || macroquad::input::is_key_down(macroquad::input::KeyCode::RightShift);
    if macroquad::input::is_key_down(macroquad::input::KeyCode::D) && shift{
        *dt += 0.001;
    }
    if macroquad::input::is_key_down(macroquad::input::KeyCode::A) && shift{
        if *dt > 0.001 {
            *dt -= 0.001;
        }
    }
}

// Increase dt when D is pressed, decrease if A is pressed
fn update_c(c: &mut f32){
    let shift = macroquad::input::is_key_down(macroquad::input::KeyCode::LeftShift) || macroquad::input::is_key_down(macroquad::input::KeyCode::RightShift);
    if macroquad::input::is_key_down(macroquad::input::KeyCode::D) && !shift {
        if *c < 1.0 {
            *c += 0.01;
            if *c > 1.0 {
                *c = 1.0;
            }
        }
    }
    if macroquad::input::is_key_down(macroquad::input::KeyCode::A) && !shift{
        if *c > 0.0 {
            *c -= 0.01;
            if *c < 0.0 {
                *c = 0.0;
            }
        }
    }
}

// Switch integration method when E is pressed (options: euler, verlet)
fn update_integration_method(integration_method: &mut IntegrationMethod){
    if macroquad::input::is_key_released(macroquad::input::KeyCode::E){
        *integration_method = match integration_method {
            IntegrationMethod::Euler => IntegrationMethod::Verlet,
            IntegrationMethod::Verlet => IntegrationMethod::Euler,
        }
    }
}

fn update_disc_vec(disc_vec: &mut Vec<Disc>, disc_idx_pair_vec: &mut Vec<(usize, usize)>, entity_r: &f32, density: &f32, charge_density: &f32){
    // Check if left-click to add new disc
    if macroquad::input::is_mouse_button_released(macroquad::input::MouseButton::Left) {
        // Calculate unique ID for new disc.
        let id = match disc_vec.iter().map(|disc| disc.id).max() {
            Some(num) => num+1,
            None => 1,
        };
        
        // New disc position.
        let pos = macroquad::input::mouse_position();

        // Create new disc and add to disc_vec
        if *entity_r < 0.0 {
            disc_vec.push(Disc::new(pos.0, pos.1, 0.0, 0.0, -*entity_r, -*density, *charge_density, id));
        } else{
            disc_vec.push(Disc::new(pos.0, pos.1, 0.0, 0.0, *entity_r, *density, *charge_density, id));
        }
        update_pair_vec(disc_idx_pair_vec, disc_vec); 
    }

    // Check if right-click to remove existing disk
    if macroquad::input::is_mouse_button_released(macroquad::input::MouseButton::Right) {
        let pos = macroquad::input::mouse_position();
        let pos_vec = vec2(pos.0, pos.1);
        disc_vec.retain(|&disc| (disc.pos.distance_squared(pos_vec) > disc.r2));
        update_pair_vec(disc_idx_pair_vec, disc_vec); 
    }

    // Check if backspace to clear all discs
    if macroquad::input::is_key_released(macroquad::input::KeyCode::Backspace) {
        disc_vec.clear();
        update_pair_vec(disc_idx_pair_vec, disc_vec); 
    }
}

fn update_pair_vec(disc_idx_pair_vec: &mut Vec<(usize, usize)>, disc_vec: &Vec<Disc>) {
    let l = disc_vec.len() as i16;

    let mut pair_vec: Vec<(usize, usize)> = Vec::new();
    for idx1 in 0..l-1 {
        for idx2 in (idx1+1)..l {
            pair_vec.push((idx1 as usize, idx2 as usize));
        }
    }
    *disc_idx_pair_vec = pair_vec;
}

// Displays the provided fps on the canvas
fn display_fps(fps: &i32){
    let mut fps_display = String::new();
    fps_display.push_str("fps: "); 
    fps_display.push_str(&fps.to_string());
    draw_text(&fps_display, 0.0, 20.0, 25.0, BLACK);
}

// Displays number of discs being simulated on the canvas
fn display_entity_cnt(entity_cnt: &i32){
    let mut entity_display = String::new();
    entity_display.push_str("entities: "); 
    entity_display.push_str(&entity_cnt.to_string());
    draw_text(&entity_display, 0.0, 40.0, 25.0, BLACK);
}

// Displays pause indicator on the canvas
fn display_pause(paused: &bool){
    if *paused {
        draw_text("Paused.    (p)", screen_width() - 160.0, 20.0, 25.0, BLACK);
    }else{
        draw_text("Playing... (p)", screen_width() - 160.0, 20.0, 25.0, BLACK);
    }
}

// Displays integration method, timestep on the canvas
fn display_integration_method(integration_method: &IntegrationMethod, dt: &f32){
    let mut label = String::new();
    match *integration_method {
        IntegrationMethod::Euler => label.push_str("euler"),
        IntegrationMethod::Verlet => label.push_str("verlet"),
    };
    label.push_str("(dt=");
    let mut dt = dt.to_string(); 
    dt.truncate(6);
    label.push_str(&dt);
    label.push_str(")");
    draw_text(&label, 0.5*screen_width()-90.0, 20.0, 25.0, BLACK);
}

// Displays size of new discs on the canvas
fn display_entity_r(entity_r: &f32){
    // Display text showing disc radius
    let mut r_display = String::new();
    r_display.push_str("disc r: "); 
    r_display.push_str(&entity_r.to_string());
    draw_text(&r_display, 0.0, 60.0, 25.0, BLACK);

    // Display bar representing disc diameter
    let diam = 2.0*entity_r;
    let mut i = 0.0;

    draw_rectangle(10.0, 0.5*screen_height(), 30.0, 1.0, BLACK);
    draw_text("+", 38.0, 0.5*screen_height()-2.5, 25.0, BLACK);
    draw_text("-", 38.0, 0.5*screen_height()+12.5, 25.0, BLACK);
 
    if diam > 0.0{
        while i < diam {
            draw_rectangle(15.0, 0.5*screen_height() - diam + i, 20.0, 1.0, BLACK);
            i += 4.0;
        }
    }else {
        while i < -diam {
            draw_rectangle(15.0, 0.5*screen_height() + i, 20.0, 1.0, BLACK);
            i += 4.0;
        }
    }
}

// Displays size of new discs on the canvas
fn display_c(c: &f32){
    // Display text showing disc radius
    let mut c_display = String::new();
    c_display.push_str("c_r: "); 
    c_display.push_str(&((*c*100.0) as i32).to_string());
    c_display.push_str("%");
    draw_text(&c_display, 0.0, screen_height()-35.0, 25.0, BLACK);

    let mut i = 0.0;
    draw_rectangle(25.0, screen_height()-29.0, 1.0, 28.0, BLACK);
    let bar_len = (screen_width()-50.0)*c;
 
    while i < bar_len {
        draw_rectangle(i+25.0, screen_height()-25.0, 1.0, 20.0, BLACK);
        i += 4.0;
    }

}

// Display help guide on the canvas
fn display_help(){
    draw_text("Click anywhere!", 0.5*screen_width()-90.0, 0.5*screen_height()-120.0, 25.0, BLACK);
    draw_text("#########################", 0.5*screen_width()-140.0, 0.5*screen_height()-100.0, 25.0, BLACK);
    draw_text("w:         increase size", 0.5*screen_width()-140.0, 0.5*screen_height()-80.0, 25.0, BLACK);
    draw_text("s:         decrease size", 0.5*screen_width()-140.0, 0.5*screen_height()-60.0, 25.0, BLACK);
    draw_text("L-clik:    add entity", 0.5*screen_width()-140.0, 0.5*screen_height()-40.0, 25.0, BLACK);
    draw_text("R-clik:    del entity", 0.5*screen_width()-140.0, 0.5*screen_height()-20.0, 25.0, BLACK);
    draw_text("backspace: clear entities", 0.5*screen_width()-140.0, 0.5*screen_height()+0.0, 25.0, BLACK);
    draw_text("e:         integr. method", 0.5*screen_width()-140.0, 0.5*screen_height()+20.0, 25.0, BLACK);
    draw_text("d:         increase c_r", 0.5*screen_width()-140.0, 0.5*screen_height()+40.0, 25.0, BLACK);
    draw_text("a:         decrease c_r", 0.5*screen_width()-140.0, 0.5*screen_height()+60.0, 25.0, BLACK);
    draw_text("p:         pause/unpause", 0.5*screen_width()-140.0, 0.5*screen_height()+80.0, 25.0, BLACK);
}

#
[macroquad::main("Gravity")]
async fn main() {    

    let mut entity_r = R_INIT;
    let mut integration_method = INTEGRATION_METHOD_INIT;

    let mut density = 1.0;
    let mut charge_density = 0.0;
    let mut c = C_INIT;
    let mut dt = DT_INIT;
    let mut fps: Vec<i32> = Vec::new();

    // vec containing all discs
    let mut disc_vec: Vec<Disc> = Vec::new();
    // vec containing disc index pairs
    let mut disc_idx_pair_vec: Vec<(usize, usize)> = Vec::new();
    
    let mut after_backstep = false;
    let mut reversal = 0.0;
    
    let mut paused: bool = false;

    loop {
        clear_background(WHITE);

        if !paused {
            let mut t_elapsed: f32 = 0.0;
            let min_dt: f32 = MAX_SLOWDOWN*dt;

            while t_elapsed < DT_PER_DISPLAY {
                let mut current_dt: f32 = dt;

                // calculate all collisions with times and retain collision events occuring before dt step
                let collisions = disc_idx_pair_vec.par_iter()
                .filter_map(|pair| { 
                    match disc_vec[pair.0].time_until_collision(&disc_vec[pair.1], &integration_method) { 
                        Some(collision_time) => if collision_time.abs() < dt {Some(((pair.0, pair.1), collision_time))} else {None}, 
                        None => None
                    }
                }).collect::<Vec<((usize, usize), f32)>>();

                // calculate soonest collision
                let soonest_collision = collisions.iter().reduce(|pair1, pair2| if pair2.1 < pair1.1 {pair2} else {pair1});

                // determine dt                    
                match soonest_collision {
                    None => {},
                    Some(mut collision) => {
                        if collision.1 < 0.0 {
                            collision = collisions.iter().filter(|pair| pair.1 < 0.0).reduce(|pair1, pair2| if pair2.1 > pair1.1 {pair2} else {pair1}).unwrap();
                        }
                        if collision.1.abs() < min_dt {
                            // if soonest_collision occurs below resolution of minimum time-step, 
                            // update timestep to minimum timestep (in the appropriate direction)
                            current_dt = collision.1.signum()*min_dt;
                        } else {
                            // otherwise set timestep to yield first collision
                            current_dt = collision.1;
                            // println!("({}",  disc_vec[collision.0.0].r+ disc_vec[collision.0.1].r - disc_vec[collision.0.0].distance(&disc_vec[collision.0.1]))
                        };
                    }
                }

                // prevent backstepping from reversing time direction
                // (limits freezing in high collision scenarios)
                if reversal+current_dt < -dt {
                    current_dt = dt;
                    reversal = 0.0;
                } else if current_dt < 0.0 {
                    reversal += current_dt;
                } else {
                    reversal = 0.0;
                }                

                match integration_method {
                    IntegrationMethod::Euler => {
                        // x(t+dt) = x(t) + v(t)*dt
                        disc_vec.par_iter_mut().for_each(|disc| disc.euler_pos_step(&current_dt));

                        // v(t+dt) = v(t) + (1/2)(F(t)/m)*dt
                        let disc_vec_cl = disc_vec.clone();
                        // note that the input here is force at x(t+dt) which does not impact velocity on this step
                        disc_vec.par_iter_mut().for_each(|disc1| disc1.euler_vel_step(&force_sum(&gravity, &disc1, &disc_vec_cl), &current_dt));
                    },
                    IntegrationMethod::Verlet => {                        
                        // update position as appropriate
                        // must be done before updating velocity bc current velocities are how collisions
                        // were calculated
                        // x(t+dt) = x(t) + v(t)*dt + (1/2)(F(t)/m)dt**2
                        disc_vec.par_iter_mut().for_each(|disc| disc.verlet_pos_step(&current_dt));

                        // apply forces & update velocity
                        // v(t+dt) = v(t) + (1/2)((F(t)+F(t+dt))/m)*dt
                        let disc_vec_cl = disc_vec.clone();
                        disc_vec.par_iter_mut().for_each(|disc1| disc1.verlet_vel_step(
                            &force_sum(&gravity, disc1, &disc_vec_cl), &current_dt));
                    }
                };

                if current_dt.abs() == min_dt {
                    // handle all collisions that occur below min dt
                    // these collisions can interact in unpredictable ways under the 
                    // resolution of the time-step, but we've decided not to care
                    collisions.iter().for_each(|pair| {
                        if pair.1 <= min_dt {
                            let (dv1, dv2) = inelastic_collision(disc_vec[pair.0.0], disc_vec[pair.0.1], c);
                            disc_vec[pair.0.0].vel += dv1;
                            disc_vec[pair.0.1].vel += dv2;
                        }
                    })
                } else if current_dt < dt {
                    // otherwise solve soonest_collision (which occurs at reasonable resolution),
                    // and dt can be updated to exactly when it happens

                    // collision handling should be PERFECT using this timestep 
                    // (assuming we don't get extremely unlucky with two simultaneous collisions)
                    let first_collision = soonest_collision.unwrap();
                    let (dv1, dv2) = inelastic_collision(disc_vec[first_collision.0.0], disc_vec[first_collision.0.1], c);
                    disc_vec[first_collision.0.0].vel += dv1;
                    disc_vec[first_collision.0.1].vel += dv2;

                    // println!("{}", disc_vec[first_collision.0.0].r+ disc_vec[first_collision.0.1].r - disc_vec[first_collision.0.0].distance(&disc_vec[first_collision.0.1]));
                }
                // println!("{}: {}", backsteps, current_dt);
                t_elapsed += current_dt;
            }
        }

        update_r(&mut entity_r);
        update_dt(&mut dt);
        update_c(&mut c);
        update_integration_method(&mut integration_method); 
        update_disc_vec(&mut disc_vec, &mut disc_idx_pair_vec, &entity_r, &density, &charge_density);

        if macroquad::input::is_key_released(macroquad::input::KeyCode::P) {
            paused = !paused;
        }

        fps.push(get_fps());
        let fps = fps.iter().rev().take(FPS_WINDOW).map(|x| *x).collect::<Vec<i32>>();
        let fps_sum: i32 = fps.iter().sum();
        let fps_cnt: i32 = fps.len() as i32;
        let fps_avg: f32 = (fps_sum as f32)/(fps_cnt as f32);

        disc_vec.iter().for_each(|disc| disc.draw());
        display_fps(&(fps_avg as i32));
        display_entity_cnt(&(disc_vec.len() as i32));
        display_entity_r(&entity_r);
        display_c(&c);
        display_pause(&paused);
        display_integration_method(&integration_method, &dt);
        if disc_vec.len() == 0 {
            display_help();
        }

        next_frame().await;
    }
}