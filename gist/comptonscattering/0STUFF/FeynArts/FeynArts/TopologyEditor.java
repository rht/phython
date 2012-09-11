/*
	TopologyEditor.java
		the topology editor UI, called by Shape
		last modified 11 Jun 03 th
*/

package de.FeynArts;

import java.math.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;
import com.wolfram.jlink.*;


public class TopologyEditor extends Frame {

  private final static double DEFAULT_RADIUS = 1.3;

  private final static int BORDER = 7, BOXWIDTH = 10, BOXHEIGHT = 10;
  private final static boolean PAINT = true, ERASE = false;

  private final static BasicStroke narrow = new BasicStroke();
  private final static BasicStroke wide = new BasicStroke(3);
  private final static BasicStroke xwide = new BasicStroke(4);

  private final static int COMMITTED = 0, CANCELLED = 1, ABORTED = 2;

  private Graphics2D screenGC;
  private Color bg;
  private double unitX, unitY;
  private int ybase;

  private CheckboxGroup positioning;
  private Checkbox gridposCheckbox;
  private Checkbox anyposCheckbox;
  private EditorCanvas editor;

  private Clickbox vertices[];
  private Propagator propagators[];
  private int nVertices, nPropagators = 0;

  private int scaleX(double x) {
    return (int)Math.round(unitX*x) + BORDER;
  }

  private int scaleY(double y) {
    return ybase - (int)Math.round(unitY*y);
  }

  private static double hypot(double x, double y) {
    return Math.sqrt(x*x + y*y);
  }

  private abstract class Clickbox {
    double x, y, backup1, backup2;
    Point center;
    Rectangle box;

    Clickbox() {
      center = new Point();
      box = new Rectangle(BOXWIDTH, BOXHEIGHT);
    }

    Clickbox(double x, double y) {
      this();
      backup1 = this.x = x;
      backup2 = this.y = y;
    }

    void update() {
      center.x = scaleX(x);
      center.y = scaleY(y);
      box.x = center.x - BOXWIDTH/2;
      box.y = center.y - BOXHEIGHT/2;
    }

    abstract void dragTo(Graphics2D g, double x, double y);
    abstract void revert();
  }

  private class Vertexbox extends Clickbox {

    Vertexbox(double x, double y) {
      super(x, y);
    }

    void revert() {
      x = backup1;
      y = backup2;
    }

    void dragTo(Graphics2D g, double x, double y) {
      double dx = x - this.x;
      double dy = y - this.y;
      if( dx == 0 && dy == 0 ) return;
      g.setColor(bg);
      g.fill(box);
      for( int i = 0; i < nPropagators; ++i ) {
        Propagator p = propagators[i];
        if( this == p.from || this == p.to ) {
          p.draw(g, ERASE);
          if( p.tadpole ) {		// a tadpole "follows" its vertex
            p.mid.x += dx;
            p.mid.y += dy;
          }
        }
      }
      this.x = x;
      this.y = y;
      update();
      for( int i = 0; i < nPropagators; ++i ) {
        Propagator p = propagators[i];
        if( this == p.from || this == p.to ) {
          p.update();
          p.draw(g, PAINT);
        }
      }
      g.setColor(Color.red);
      g.fill(box);
    }
  }

  private class Propagatorbox extends Clickbox {
    Propagator parent;
    private double backupheight;

    Propagatorbox(double x, double y, Propagator parent) {
      super(x, y);
      this.parent = parent;
      backupheight = parent.height;
    }

    void revert() {
      x = backup1;
      y = backup2;
      parent.height = backupheight;
      parent.update();
    }

    void dragTo(Graphics2D g, double x, double y) {
      if( x == this.x && y == this.y ) return;
      if( parent.tadpole ) parent.draw(g, ERASE);
      else {
        double xm = .5*(parent.to.x + parent.from.x);
        double ym = .5*(parent.to.y + parent.from.y);
        double dx = parent.to.x - parent.from.x;
        double dy = parent.to.y - parent.from.y;

	// adjust clicked position to lie on the perpendicular bisector
        if( dx == 0 ) y = ym;
        else if( dy == 0 ) x = xm;
        else {
          double h = dy/dx;
          double d = h + dx/dy;
          x = (xm/h + ym - y + h*x)/d;
          y = (xm + ym*h - x + y/h)/d;
        }

        double height = hypot(xm - x, ym - y);
        if( height < .3 ) height = 0;
        else {
          height = 2*height/hypot(dx, dy);
	    // cross product tells which side of the prop we are on
          double cross = dy*(x - parent.from.x) - dx*(y - parent.from.y);
          if( height*cross < 0 ) height = -height;
        }
        if( height == parent.height ) return;
        parent.draw(g, ERASE);
        if( (height >= 0 && parent.height < 0) ||
            (height < 0 && parent.height >= 0) )
          parent.label.angle = -parent.label.angle;
        parent.height = height;
      }
      this.x = x;
      this.y = y;
      parent.update();
      parent.draw(g, PAINT);
    }
  }

  private class Labelbox extends Clickbox {
    Propagator parent;
    double radius, angle;

    Labelbox(double radius, double angle, Propagator parent) {
      super();
      this.parent = parent;
      backup1 = this.radius = radius;
      backup2 = this.angle = angle;
    }

    void revert() {
      radius = backup1;
      angle = backup2;
    }

    void setXY() {
      double phi = angle + parent.arcCenterAngle;
      x = parent.mid.x + radius*Math.cos(phi);
      y = parent.mid.y + radius*Math.sin(phi);
    }

    void dragTo(Graphics2D g, double x, double y) {
      if( x == this.x && y == this.y ) return;
      double dx = x - parent.mid.x;
      double dy = y - parent.mid.y;
      radius = hypot(dx, dy);
      angle = Math.atan2(dy, dx) - parent.arcCenterAngle;
      if( Math.abs(angle) < .1 ) angle = 0;
      else if( Math.abs(angle - Math.PI) < .1 ) angle = Math.PI;

      g.setColor(bg);
      g.drawLine(parent.mid.center.x, parent.mid.center.y, center.x, center.y);
      g.fill(box);
      setXY();
      update();
      g.setColor(Color.black);
      g.drawLine(parent.mid.center.x, parent.mid.center.y, center.x, center.y);
      g.setColor(Color.green);
      g.fill(box);
    }
  }


  private class Propagator {
    Vertexbox from, to;
    Propagatorbox mid;
    Labelbox label;
    boolean tadpole;
    double height, arcCenterAngle;
    private double arcRadius, arcCenterX, arcCenterY;
    private int arcAngleStart, arcAngleExtent;

    Propagator(double prop[], int off) {
      if( tadpole = prop[off] < 0 ) {
        from = to = (Vertexbox)vertices[-(int)prop[off] - 1];
        mid = new Propagatorbox(prop[off + 1], prop[off + 2], this);
      }
      else {
        from = (Vertexbox)vertices[(int)prop[off] - 1];
        to = (Vertexbox)vertices[(int)prop[off + 1] - 1];
        height = prop[off + 2];
        mid = new Propagatorbox(0, 0, this);
      }
      label = new Labelbox(prop[off + 3], prop[off + 4], this);
      update();
    }

    void update() {
      if( tadpole ) {
        arcCenterX = .5*(from.x + mid.x);
        arcCenterY = .5*(from.y + mid.y);
        double dx = arcCenterX - from.x;
        double dy = arcCenterY - from.y;
        arcRadius = hypot(dx, dy);
        arcCenterAngle = Math.atan2(dy, dx);
        arcAngleStart = 0;
        arcAngleExtent = 360;
      }
      else {
        double dx = to.x - from.x;
        double dy = to.y - from.y;
        double halfLength = .5*hypot(dx, dy);
        double centerAngle = Math.atan2(dy, dx) - .5*Math.PI;
        double halfAngle, rad;
        arcCenterX = mid.x = .5*(from.x + to.x);
        arcCenterY = mid.y = .5*(from.y + to.y);
        if( height == 0 ) {
          arcRadius = rad = 20000;
          halfAngle = Math.asin(halfLength/rad);
        }
        else {
          if( height < 0 ) centerAngle += Math.PI;
          halfAngle = 2*Math.atan(rad = Math.abs(height));
          mid.x += rad*halfLength*Math.cos(centerAngle);
          mid.y += rad*halfLength*Math.sin(centerAngle);
          arcRadius = rad = halfLength/Math.sin(halfAngle);
          rad = Math.sqrt(rad*rad - halfLength*halfLength);
          if( Math.abs(height) > 1 ) rad = -rad;
        }
        arcCenterX -= rad*Math.cos(centerAngle);
        arcCenterY -= rad*Math.sin(centerAngle);
        arcCenterAngle = centerAngle;
        arcAngleStart = (int)Math.round(Math.toDegrees(centerAngle - halfAngle));
        arcAngleExtent = (int)Math.round(Math.toDegrees(2*halfAngle));
      }
      label.setXY();
    }

    void draw(Graphics2D g, boolean inColor) {
      if( inColor ) {
        g.setStroke(wide);
        g.setColor(Color.black);
      }
      else {
        g.setStroke(xwide);
        g.setColor(bg);
      }
      if( arcRadius < 20000 )
        g.drawArc(
          scaleX(arcCenterX - arcRadius),
          scaleY(arcCenterY + arcRadius),
          (int)Math.round(2*arcRadius*unitX),
          (int)Math.round(2*arcRadius*unitY),
          arcAngleStart, arcAngleExtent);
      else
        g.drawLine(from.center.x, from.center.y, to.center.x, to.center.y);
      g.setStroke(narrow);
      mid.update();
      label.update();
      g.drawLine(mid.center.x, mid.center.y, label.center.x, label.center.y);
      if( inColor ) g.setColor(Color.blue);
      g.fill(mid.box);
      if( inColor ) g.setColor(Color.green);
      g.fill(label.box);
    }
  }


  public class EditorCanvas extends Canvas
    implements ComponentListener, MouseListener, MouseMotionListener {

    private Clickbox selected;
    private Image current = null;
    private Graphics2D currentGC;
    private boolean redraw = true;

    EditorCanvas(int width, int height) {
      setSize(width, height);
      ybase = height - BORDER;
      unitX = (width - 2*BORDER)/20.;
      unitY = (height - 2*BORDER)/20.;
      addComponentListener(this);
      addMouseListener(this);
      addMouseMotionListener(this);
    }

    public void update(Graphics g) {
      paint(g);
    }

    public void paint(Graphics g) {
      if( nPropagators == 0 ) return;
      if( current == null || redraw ) {
        Dimension size = getSize();
        if( current == null ) {
          current = createImage(size.width, size.height);
          currentGC = (Graphics2D)(current == null ? g : current.getGraphics());
        }
        currentGC.setColor(getBackground());
        currentGC.fillRect(0, 0, size.width, size.height);
        paintCanvas(currentGC);
        redraw = false;
      }
      if( current != null ) g.drawImage(current, 0, 0, null);
    }

    private void paintCanvas(Graphics2D g) {
      g.setColor(Color.black);
      double x = BORDER;
      for( int i = 20; i >= 0; --i ) {		// draw grid
        int i5 = i % 5;
        double y = ybase;
        for( int j = 20; j >= 0; --j ) {
          int ix = (int)Math.round(x);
          int iy = (int)Math.round(y);
          if( i5 == 0 && j % 5 == 0 ) {
            g.drawLine(ix - 2, iy, ix + 2, iy);
            g.drawLine(ix, iy - 2, ix, iy + 2);
          }
          else g.drawRect(ix, iy, 0, 0);
          y -= unitY;
        }
        x += unitX;
      }

      for( int i = 0; i < nVertices; ++i ) vertices[i].update();

      for( int i = 0; i < nPropagators; ++i )
        propagators[i].draw(g, PAINT);

      g.setColor(Color.red);
      for( int i = 0; i < nVertices; ++i )
        g.fill(vertices[i].box);
    }

    public void componentResized(ComponentEvent e) {
      Dimension size = getSize();
      ybase = (int)Math.round(size.height) - BORDER;
      unitX = (size.width - 2*BORDER)/20.;
      unitY = (size.height - 2*BORDER)/20.;
      currentGC.dispose();
      current = null;
    }

    public void componentMoved(ComponentEvent e) {}
    public void componentHidden(ComponentEvent e) {}
    public void componentShown(ComponentEvent e) {}

    public void mouseDragged(MouseEvent m) {
      if( selected != null ) {
        double x = (m.getX() - BORDER)/unitX;
        double y = (ybase - m.getY())/unitY;
        if( positioning.getSelectedCheckbox() == gridposCheckbox ) {
          x = .5*Math.round(2*x);
          y = .5*Math.round(2*y);
        }
        x = Math.max(Math.min(x, 20), 0);
        y = Math.max(Math.min(y, 20), 0);
        selected.dragTo(screenGC, x, y);
      }
    }

    public void mousePressed(MouseEvent m) {
      if( SwingUtilities.isLeftMouseButton(m) ) {
        for( int i = 0; i < vertices.length; ++i )
          if( vertices[i].box.contains(m.getPoint()) ) {
            selected = vertices[i];
            break;
          }
      }
    }

    public void mouseReleased(MouseEvent m) {
      if( selected != null ) {
        selected = null;
        redraw = true;
        repaint();
      }
    }

    public void mouseClicked(MouseEvent m) {
      if( SwingUtilities.isLeftMouseButton(m) ) return;

      Clickbox sel = null;
      for( int i = nVertices; i < vertices.length; ++i )
        if( vertices[i].box.contains(m.getPoint()) ) {
          sel = vertices[i];
          break;
        }
      if( sel == null ) return;

      if( sel instanceof Propagatorbox ) {
        Propagator p = ((Propagatorbox)sel).parent;
        if( p.tadpole ) return;
        if( SwingUtilities.isRightMouseButton(m) ) {
          p.height = -p.height;
          p.label.angle = -p.label.angle;
        }
        else {
          if( p.height < 0 ) p.label.angle = -p.label.angle;
          p.height = 0;
        }
        p.update();
      }
      else {
        Labelbox l = (Labelbox)sel;
        if( SwingUtilities.isRightMouseButton(m) )
          l.angle = Math.PI - l.angle;
        else {
          l.radius = DEFAULT_RADIUS;
          l.angle = 0;
        }
        l.setXY();
      }
      redraw = true;
      repaint();
    }

    public void mouseMoved(MouseEvent m) {}
    public void mouseEntered(MouseEvent m) {}
    public void mouseExited(MouseEvent m) {}
  }

  public class HelpCanvas extends Canvas {
    private final int WIDTH = 107, HEIGHT = 175;
    private Image mousehelp;

    HelpCanvas() {
      setSize(WIDTH, HEIGHT);
    }

    public void paint(Graphics g) {
      if( mousehelp == null ) {
        mousehelp = createImage(WIDTH, HEIGHT);
        Graphics gi = mousehelp.getGraphics();
        gi.setColor(getBackground());
        gi.fillRect(0, 0, WIDTH, HEIGHT);

        final int TEXTHEIGHT = 16;
        final int SPACER = 10;
        final int LINE0 = TEXTHEIGHT;
        final int LINE1 = LINE0 + SPACER + TEXTHEIGHT;
        final int LINE2 = LINE1 + TEXTHEIGHT;
        final int LINE3 = LINE2 + SPACER + TEXTHEIGHT;
        final int LINE4 = LINE3 + TEXTHEIGHT;
        final int LINE5 = LINE4 + TEXTHEIGHT;
        final int LINE6 = LINE5 + SPACER + TEXTHEIGHT;
        final int LINE7 = LINE6 + TEXTHEIGHT;
        final int INDENT0 = BORDER;
        final int INDENT1 = INDENT0 + 3 + BOXWIDTH;
        final int INDENT2 = INDENT1 + 3 + BOXWIDTH;
        final int INDENT3 = INDENT2 + 3 + BOXWIDTH;

        gi.setColor(Color.black);
        gi.drawString("mouse:", INDENT0, LINE0);
        gi.drawString("left button", INDENT0, LINE1);
        gi.drawString("drag", INDENT3, LINE2);

        gi.drawString("middle button", INDENT0, LINE3);
        gi.drawString("straighten", INDENT1, LINE4);
        gi.drawString("default pos.", INDENT1, LINE5);

        gi.drawString("right button", INDENT0, LINE6);
        gi.drawString("mirror", INDENT2, LINE7);

        gi.setColor(Color.red);
        gi.fillRect(INDENT0, LINE2 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);

        gi.setColor(Color.blue);
        gi.fillRect(INDENT1, LINE2 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);
        gi.fillRect(INDENT0, LINE4 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);
        gi.fillRect(INDENT0, LINE7 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);

        gi.setColor(Color.green);
        gi.fillRect(INDENT2, LINE2 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);
        gi.fillRect(INDENT0, LINE5 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);
        gi.fillRect(INDENT1, LINE7 - BOXHEIGHT + 1, BOXWIDTH, BOXHEIGHT);
      }

      g.drawImage(mousehelp, 0, 0, null);
    }
  }

  public TopologyEditor() {
    super("FeynArts Topology Editor");

    editor = new EditorCanvas(20*20 + 2*BORDER, 20*20 + 2*BORDER);

    Box rightpanel = Box.createVerticalBox();

    Button okButton = new Button("OK");
    okButton.addActionListener(
      new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          quit(COMMITTED);
        }
      });
    rightpanel.add(okButton);

    Button cancelButton = new Button("Cancel");
    cancelButton.addActionListener(
      new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          quit(CANCELLED);
        }
      });
    rightpanel.add(cancelButton);

    Button abortButton = new Button("Abort");
    abortButton.addActionListener(
      new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          quit(ABORTED);
        }
      });
    rightpanel.add(abortButton);

    Button revertButton = new Button("Revert");
    revertButton.addActionListener(
      new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          for( int i = 0; i < vertices.length; ++i )
            vertices[i].revert();
          editor.redraw = true;
          editor.repaint();
        }
      });
    rightpanel.add(revertButton);

    rightpanel.add(Box.createVerticalStrut(10));

    positioning = new CheckboxGroup();
    gridposCheckbox = new Checkbox("Grid position", positioning, true);
    anyposCheckbox = new Checkbox("Any position", positioning, false);
    rightpanel.add(gridposCheckbox);
    rightpanel.add(anyposCheckbox);

    rightpanel.add(Box.createVerticalStrut(10));

    HelpCanvas help = new HelpCanvas();
    rightpanel.add(help);

    add(BorderLayout.CENTER, editor);
    add(BorderLayout.EAST, rightpanel);
    pack();
    addWindowListener(new WindowAdapter() {
      public void windowClosing(WindowEvent e) {
        quit(ABORTED);
      }
    });
    setVisible(true);

    screenGC = (Graphics2D)editor.getGraphics();
    bg = editor.getBackground();
  }

  public void closeWindow() {
    dispose();
    screenGC.dispose();
  }

  private void quit(int success) {
    KernelLink link = StdLink.getLink();
    if( link == null ) return;
    try {
      link.evaluate("(EndModal[]; " + success + ")");
      link.discardAnswer();
    }
    catch( MathLinkException e ) {
      System.err.println(e.getMessage());
    }
  }

  public void putShapeData(double vert[], double prop[]) {
    vertices = new Clickbox[vert.length/2 + 2*prop.length/5];
    propagators = new Propagator[prop.length/5];

    int n = 0;
    for( int i = 0; i < vert.length; i += 2 )
      vertices[n++] = new Vertexbox(vert[i], vert[i + 1]);
    nVertices = n;

    nPropagators = 0;
    for( int i = 0; i < prop.length; i += 5 ) {
      Propagator p = new Propagator(prop, i);
      propagators[nPropagators++] = p;
      vertices[n++] = p.label;
      vertices[n++] = p.mid;
	// note: the order "first label, then mid" is important
	//       because mid.revert() calls label.update()
    }
    editor.redraw = true;
    editor.repaint();
  }

  public void getShapeData() {
    KernelLink link = StdLink.getLink();
    if( link == null ) return;

    link.beginManual();
    try {
      link.putFunction("List", 3);

      link.putFunction("List", nVertices);
      for( int i = 0; i < nVertices; ++i ) {
        link.putFunction("List", 2);
        link.put(vertices[i].x);
        link.put(vertices[i].y);
      }

      link.putFunction("List", nPropagators);
      for( int i = 0; i < nPropagators; ++i ) {
        Propagator p = propagators[i];
        if( p.tadpole ) {
          link.putFunction("List", 2);
          link.put(p.mid.x);
          link.put(p.mid.y);
        }
        else link.put(p.height);
      }

      link.putFunction("List", nPropagators);
      for( int i = 0; i < nPropagators; ++i ) {
        Labelbox l = propagators[i].label;
        if( l.angle == 0 )
          link.put(l.radius/DEFAULT_RADIUS);
        else if( l.angle == Math.PI )
          link.put(-l.radius/DEFAULT_RADIUS);
        else {
          link.putFunction("List", 2);
          link.put(l.radius);
          link.put(l.angle);
        }
      }
//      link.endPacket();
    }
    catch( MathLinkException e ) {
      System.err.println(e.getMessage());
    }
  }

}
